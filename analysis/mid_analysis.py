#%matplotlib inline
import theano
theano.config.floatX = 'float64'
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pymc3 as pm
import pandas as pd
from glob import glob
import ipdb

def get_subj_num(filename):
    f = open(filename, 'r')
    header = f.readline()
    f.close()
    subj_num = int(header[header.find('ID:')+4:header.find('Task type')])
    return subj_num

def prep_data(filename):
    subj_num = get_subj_num(filename)
    df = pd.read_csv(filename, delimiter='\t', header=1)
    df = df[~pd.isnull(df.Trial)]
    
    df.columns = df.columns.str.strip()
    df = df.drop(df.columns[df.columns.str.find('Unna')>-1], 1) # drop extraneous columns

    # calculate RT and code by side (left: 0, right: 1)
    df['response'] = (df['Response Made by Subject'].str.find('Right') > -1).astype(int)
    df['rt'] = (df.loc[:, 'Response time'] - df.loc[:, 'Target Phase Start Time']) * 10e-4
    df = df[~pd.isnull(df.rt)] # remove trials with no response

    # get the pavlovian bias terms (need prev correct as well as nuisence regressor)
    df['prev_correct'] = np.concatenate([[0], np.array(df.Outcome == 'SUCCESS', dtype=int)[:-1]])
    df['prev_reward'] = np.concatenate([[0], np.array(df.Amount>0, dtype=int)[:-1]])

    df['subj_idx'] = [subj_num] * len(df)
    df['trial_type'] = df['Trial Category']
    df['target_duration'] = df['Target Phase Duration'] * 10e-4
    
    # Code the target direction (as well as possible, this is a guess)
    df['stim'] = ['N/A'] * len(df)
    R = (df['Response Made by Subject'].str.find('Left') == 0) & (df['Outcome']=="Failure") | \
        (df['Response Made by Subject'].str.find('Right') == 0 & (df['Outcome']=="Success")) | \
        (df['Response Made by Subject'].str.find('Right') > 0)
    L = (df['Response Made by Subject'].str.find('Right') == 0) & (df['Outcome']=="Failure") | \
        (df['Response Made by Subject'].str.find('Left') == 0 & (df['Outcome']=="Success")) | \
        (df['Response Made by Subject'].str.find('Left') > 0)
    df.loc[R, 'stim'] = 'R'
    df.loc[L, 'stim'] = 'L'
    df['subj_idx'] = [subj_num] * len(df)
#     if sum(df['rt'] > 0) == 0:
#         print df
#         raise
    df = df[df['rt'] > 0] # restrict sample to valid trials
    
    return df[
        ['stim', 'rt', 'trial_type', 'target_duration', 
         'prev_correct', 'prev_reward', 'subj_idx']
    ]

### Read files
files = glob('/home/dan/projects/imagen/data/archived/Recieved/MID/preprocessed/MID_BL_behavioraldata/mid_*.csv')
subj_list = pd.read_csv('~/projects/imagen/data/sandbox_subject_list.csv').Subject.values
subj_list.sort()

all_data = []
for ii, f in enumerate(files):
    if get_subj_num(f) in subj_list:
      print('Reading file: ' + f)
      all_data.append(prep_data(f))

# Task conditions - associated w/ indicator variables
conditions = ('disp_side', 'tau_target', 'reward', 'high_reward', \
   'low_reward', 'prev_reward', 'prev_correct')

# Data is comprised of 3 book-keeping variables and indicators for each condition
group_data = {el:[] for el in conditions}
group_data['rt'] = []
group_data['subj_id'] = []
group_data['trial_category'] = []

t = 0
for ii, data in enumerate(all_data):
    if data.shape[0] > 0:
        group_data['subj_id']        += [t] * data.shape[0]
        #group_data['rt']             += list(1 / np.array(data['rt']))
        group_data['rt']             += list(np.array(data['rt']))
        group_data['disp_side']      += np.array(data.stim.str.find('L') == 0, dtype=int).tolist()
        group_data['tau_target']     += list(1 / np.array(data['target_duration']) - np.mean(1 / np.array(data['target_duration'])))
        group_data['reward']         += np.array(data.trial_type.str.find('NO_WIN') == 0, dtype=int).tolist()
        group_data['low_reward']     += np.array(data.trial_type.str.find('SMALL_WIN') > -1, dtype=int).tolist()
        group_data['high_reward']    += np.array(data.trial_type.str.find('BIG_WIN') > -1, dtype=int).tolist()
        group_data['prev_correct']   += data.prev_correct.tolist()
        group_data['prev_reward']    += data.prev_reward.tolist()
        group_data['trial_category'] += data.trial_type.tolist()
        t += 1

group_data = pd.DataFrame(group_data)
n_subj = len(set(group_data.subj_id))

print("Subject's Requested: %d" % len(subj_list))
print("Subject's In dataset: %d" % n_subj)

### Display cumulative distribtions of 1/RT
#fig, ax = plt.subplots(figsize=(8, 6))
#cc = sns.color_palette("Set2")
#summary = list()
#for label, color in zip('BIG_WIN SMALL_WIN NO_WIN'.split(), cc):
#    
#    x = group_data.loc[group_data.trial_category==label, 'rt']
#    x = x[~np.isnan(x)]
#    x.sort()
#    x = np.flipud(x)
#    rt = np.cumsum(np.ones(len(x)) / len(x))
#    plt.plot(x, rt, 'd', color=color, label=label)
#
#ax.set_xlim(8, 1)
#plt.legend(loc='lower right')
#ax.set_yscale('logit')
#ax.set_xlabel('RT (Hz)')
#plt.show()

### Define Hierarchical Models ###

# Group model: subject ID and each condition
subj_idx     = group_data['subj_id'].values   
tau_target   = group_data.tau_target.values
disp_side    = group_data.disp_side.values
reward       = group_data.reward.values
low_reward   = group_data.low_reward.values
high_reward  = group_data.high_reward.values
prev_correct = group_data.prev_correct.values
prev_reward  = group_data.prev_reward.values

# Plausible models

### Ex-Gaussian Model ###
model_conds = conditions
with pm.Model() as later_model:

   # Using dicts for flexible condition inclusion
   group_params = dict()
   subj_params  = dict()

   # Hyperpriors for group nodes
   group_params['mu_intercept']  = pm.Normal('mu_intercept', mu = 450, sd = 100)
   #group_params['sd_intercept'] = pm.Normal('sd_intercept', mu = 0, sd = 100**2)
   group_params['sd_intercept']  = pm.Uniform('sd_intercept', lower = 0, upper = 200)

   for cond in model_conds:
      group_params['mu_beta_'+cond] = pm.Normal( 'mu_beta_' + cond, mu = 0., sd = 100**2)
      group_params['sd_beta_'+cond] = pm.Uniform('sd_beta_' + cond, lower = 0, upper = 200)
      #group_params['sd_beta_'+cond] = pm.HalfCauchy('sd_beta_' + cond, beta = 100)

   # Hyperprior for each subjects' scale parameter
   scale_stdev = pm.HalfCauchy('scale_stdev', beta=100)


   # Subject level parameters
   subj_params['intercept'] = pm.Normal('intercept', mu = group_params['mu_intercept'], sd = group_params['sd_intercept'], shape = n_subj)

   for cond in model_conds:
      name = 'beta_' + cond
      subj_params[name] = pm.Normal(name, mu = group_params['mu_'+name], sd = group_params['sd_'+name], shape = n_subj)


   # Model error
   #stdev     = pm.Uniform('eps', lower=0, upper=100, shape=n_subj)
   stdev     = pm.HalfCauchy('eps', beta = 100, shape = n_subj)
   stdev_all = stdev[subj_idx] 

   nu     = pm.Uniform('tau', lower = 0, upper = 50, shape = n_subj)
   nu_all = nu[subj_idx] 

   # Regression Model
   y_hat = subj_params['intercept'][subj_idx]
   #for cond in model_conds: 
   #   y_hat = y_hat + subj_params['beta_'+cond][subj_idx]*group_data[cond].values
    
   # Data likelihood
   y_like = pm.ExGaussian('like', mu = y_hat, sigma = stdev_all, nu = nu_all, observed = group_data.rt.values)


#print('Unimputed censored model')
#with unimputed_censored_model:
#  trace = pm.sample(n)
#  pm.plot_posterior(trace[-1000:], varnames=['mu', 'sigma'])
#  plt.show()

niter = 100
njobs = 1

### Model fitting ###
def fit_model(model):

    #with model:
    #    means, sds, elbos = pm.variational.advi(n = 100)

    # Inference button (TM)!
    with model:
        start = pm.find_MAP()
        step  = pm.NUTS()
        trace = pm.sample(niter, step, start = start, njobs = njobs)
        
    return trace

trace_1 = fit_model(later_model)





### Model comparisons ###
#traces = [trace_1]
#traces = [hierarchical_trace_1, hierarchical_trace_2, hierarchical_trace_3, 
#          hierarchical_trace_4, hierarchical_trace_5, hierarchical_trace_6,
#          hierarchical_trace_7, hierarchical_trace_8, trace_1,
#          hierarchical_trace_10,
#         ]

#models = [later_model]
#models = [hierarchical_later_model, hierarchical_model_2, hierarchical_model_3, 
#          hierarchical_model_4, hierarchical_model_5, hierarchical_model_6,
#          hierarchical_model_7, hierarchical_model_8, hierarchical_model_9,
#          hierarchical_later_model0,
#         ]
#
#fit_statitics = pd.DataFrame(index=['M%d' % ii for ii in range(1, len(models)+1)], 
#                             columns=['DIC', 'WAIC'])
#fit_statitics.index.name = 'Model'
#
#for ii, (trace, model) in enumerate(zip(traces, models), 1):
#    fit_statitics.loc["M%d" % ii, 'DIC' ] = pm.stats.dic( model = model, trace = trace[500:])
#    fit_statitics.loc["M%d" % ii, 'WAIC'] = pm.stats.waic(model = model, trace = trace[500:])
#
#
#
#_fit_statitics = pd.melt(fit_statitics.reset_index(), id_vars=['Model'], 
#                         var_name='fit_measure', value_name='statistic')
#
#g = sns.factorplot(x='Model', rt='statistic', col='fit_measure', hue='fit_measure', 
#                   data=_fit_statitics, 
#                   kind='bar')
#
#_fit_statitics
#
#
#tmp = _fit_statitics[_fit_statitics.fit_measure=='WAIC']
#tmp[tmp.statistic==tmp.statistic.min()]
#
#
#trace_1.varnames
#
pm.traceplot(trace_1, varnames=['intercept'])
plt.show()

pm.autocorrplot(trace_1, varnames=['intercept'])
plt.show()


pm.effective_n(trace_1)
pm.gelman_rubin(trace_1)

plt.plot(pm.geweke(trace_1['intercept'])[:,1], 'o')
plt.axhline(1, c='red')
plt.axhline(-1, c='red')
plt.gca().margins(0.05)

pm.summary(trace_1, varnames=['intercept'])

with coin_context:
    ppc = pm.sample_ppc(trace_1, samples=100)
ppc['y'].ravel()

#niter = 2000
#with pm.Model() as text_save_demo:
#    p = pm.Beta('p', alpha=2, beta=2)
#    y = pm.Binomial('y', n=n, p=p, observed=heads)
#    db = Text('trace')
#    trace = pm.sample(niter, trace=db)



pgm = daft.PGM(shape=[2.5, 3.0], origin=[0, -0.5])

pgm.add_node(daft.Node("alpha", r"$\alpha$", 0.5, 2, fixed=True))
pgm.add_node(daft.Node("beta", r"$\beta$", 1.5, 2, fixed=True))
pgm.add_node(daft.Node("p", r"$p$", 1, 1))
pgm.add_node(daft.Node("n", r"$n$", 2, 0, fixed=True))
pgm.add_node(daft.Node("y", r"$y$", 1, 0, observed=True))

pgm.add_edge("alpha", "p")
pgm.add_edge("beta", "p")
pgm.add_edge("n", "y")
pgm.add_edge("p", "y")

pgm.render()







from pymc3.backends import SQLite

niter = 2000
with pm.Model() as sqlie3_save_demo:
    p = pm.Beta('p', alpha=2, beta=2)
    y = pm.Binomial('y', n=n, p=p, observed=heads)
    db = SQLite('trace.db')
    trace = pm.sample(niter, trace=db)

with sqlie3_save_demo:
    trace = pm.backends.sqlite.load('trace.db')
    pm.traceplot(trace, varnames=['p'])


with pm.Model() as prior_context:
    sigma = pm.Gamma('gamma', alpha=2.0, beta=1.0)
    mu = pm.Normal('mu', mu=0, sd=sigma)
    trace = pm.sample(niter, step=pm.Metropolis())

sns.distplot(xs, rug=True)


####
# Snippets
####
# Unimputed censored model
# Helper function: Normal left-censored CDF
#def normal_lcdf(mu, sigma, x):
#    z = (x - mu) / sigma
#    return tt.switch(
#        tt.lt(z, -1.0),
#        tt.log(tt.erfcx(-z / tt.sqrt(2.)) / 2.) - tt.sqr(z) / 2.,
#        tt.log1p(-tt.erfc(z / tt.sqrt(2.)) / 2.)
#    )
#
#def censored_left_likelihood(mu, sigma, n_left_censored, lower_bound):
#    return n_left_censored * normal_lcdf(mu, sigma, lower_bound)
#
#with pm.Model() as unimputed_censored_model:
#    mu       = pm.Normal('mu', mu=0., sd=(high - low) / 2.)
#    sigma    = pm.HalfNormal('sigma', sd=(high - low) / 2.)
#    observed = pm.Normal('observed', mu=mu, sd=sigma, observed=truncated, shape=n_subj)
#    
#    left_censored = pm.Potential( 'left_censored', censored_left_likelihood(mu, sigma, n_left_censored, low))



# Define models:
#model_conds  = conditions
#with pm.Model() as later_model:
#
#   # Using dicts for flexible condition inclusion
#   group_params = dict()
#   subj_params  = dict()
#
#   # Hyperpriors for group nodes
#   group_params['mu_intercept'] = pm.Normal('mu_intercept', mu = 0, sd = 100**2)
#   group_params['sd_intercept'] = pm.Normal('sd_intercept', mu = 0, sd = 100**2)

#   for cond in model_conds:
#      group_params['mu_beta_'+cond] = pm.Normal( 'mu_beta_' + cond, mu = 0., sd = 100**2)
#      #group_params['sd_beta_'+cond] = pm.Uniform('sd_beta_' + cond, lower = 0, upper = 100)
#      group_params['sd_beta_'+cond] = pm.HalfCauchy('sd_beta_' + cond, beta = 100)
#
#   # Hyperprior for each subjects' scale parameter
#   scale_stdev = pm.HalfCauchy('scale_stdev', beta=100)
#
#
#   # Subject level parameters
#   subj_params['intercept'] = pm.Normal('intercept', mu = group_params['mu_intercept'], sd = group_params['sd_intercept'], shape = n_subj)
#
#   for cond in model_conds:
#      name = 'beta_' + cond
#      subj_params[name] = pm.Normal(name, mu = group_params['mu_'+name], sd = group_params['sd_'+name], shape = n_subj)
#
#
#   # Model error
#   #stdev     = pm.Uniform('eps', lower=0, upper=100, shape=n_subj)
#   stdev     = pm.HalfCauchy('eps', beta = 100, shape = n_subj)
#   stdev_all = stdev[subj_idx] 
#
#   # Regression Model
#   y_hat = subj_params['intercept'][subj_idx]
#   for cond in model_conds: 
#      y_hat = y_hat + subj_params['beta_'+cond][subj_idx]*group_data[cond].values
#    
#   # Data likelihood
#   y_like = pm.Normal('like', mu = y_hat, sd = stdev_all, observed = group_data.rt.values)

