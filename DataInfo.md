# IMAGEN Misc. Info #

## Table of Contents: ##
* [Outline](#outline)
* [Summary Stats for MID and SST Params](#summary stats for mid and sst params)
* [Mean Shift and Normalize](#mean shift and normalize)
* [Shuffle Test](#Shuffle Test)
* [Further Data Subdivision](#further data subdivision)
* [Missing Data: Age 14](#missing data age 14)
* [Impute Features](#impute features)
* [Clean Up Features](#clean up features)
* [Misc. Correlations](#misc correlations)
* [Canonical Correlations Analysis: Age 14](#canonical correlations analysis age 14)
* [Principle Components Analysis](#principle components analysis)
* [CCA Across Baseline and Follow Up](#cca across baseline and follow up)
* [CCA from Age 14 Tasks to Age 18 Survey](#cca from age 14 tasks to age 18 survey)
* [Gaussian Mixture Model](#gaussian mixture model)

## Outline <a class="anchor" id="outline"></a>
**The following should be determined in this file:**
- What data exists, what subset is/ shoudl be used, and its raw state.
- How much needs to be interpolated or removed, how, and what other QC needs to be done.
    - How important statistics of the data are changed by these processes.
- What analyses make the most sense, and what we roughly expect from them

**It also includes...**
- Attempted reproductions of others' work
- Sanity checks on work done here and on the data

**Actual steps in achieving these goals include:**
- Stably inducting and partitioning the data
- 

# Misc. Info #

## Description of the CANTAB Measures ##
Participants completed five of the CANTAB tests. 

**The Affective Go/No-go** task comprised of alternating blocks in which participants were presented with positively or negatively valenced target words embedded in a stream of neutral distracter words. Participants were instructed to respond to targets with a button press.  Measures included in the analyses were the total number of omissions to positive and negative targets, and the average response latency to positive and negative target words.
- 'agn_mean_correct_latency_negative'  Measure is average RT to negative targets
- 'agn_mean_correct_latency_positive' Measure is average RT to positive targets
- 'agn_total_omissions_negative' Measure is number of negatively valenced target words that the participant failed to respond to 
- 'agn_total_omissions_positive' Measure is number of positively valenced target words that the participant failed to respond to 

<br />
In the **Pattern Recognition Memory** task participants were required to remember 12 abstract patterns; the percentage of patterns correctly recognized on a two alternative forced choice task completed immediately after encoding was included in the analyses. 
- 'prm_percent_correct'  Score reflects % correct

<br />
**The Spatial Working Memory Task** required participants to “search” for a token hidden by one of a number of boxes on the monitor by selecting the boxes in sequence. Once the token is uncovered, participants must search again with the condition that the token will not be hidden in the same location more than once. The number of times participants returned to search a box that had already contained the token was entered into the analyses as an error measure. We also included a strategy score (ranging from 1-37, with lower scores indicating a more strategic approach), which reflects how often a search sequence was initiated from a novel position.
- 'swm_between_errors':  Score reflects the number of times the child revisits a box that has already been shown to contain a token (higher scores represent less efficiency).
- 'swm_strategy': For displays with more than 6 boxes, this measure represents the number of boxes in which the child initiates a new search within the same display (higher scores represent a less strategic approach).

<br />
**The Rapid Visual Information Processing** task comprised of a stream of digits presented at 1.67Hz and participants were required to monitor the stream for target sequence of three digits. We included a signal detection measure of sensitivity to the target sequence in the analyses. 
- 'rvp_a'  Score is a signal detection measure reflecting ability to discriminate a target sequence of 3 digits embedded in a rapidly presented stream.

<br />
**The Cambridge Guessing Task (CGT)** was a modified version of the Cambridge Gambling Task, renamed in order to make it appropriate to administer to adolescents.  On each trial of the CGT the participant was presented with 10 boxes, some of which are blue, some of which are red, and must “guess” which color box conceals a hidden yellow token.  Participants start the task with 100 points and lose or acquire points by wagering on their guess.  The options the participant can choose to wager are determined by the program as a proportion of their total number of points, presented in either increasing or decreasing amounts.  Measures included the time taken to select the option on which to bet, an average of the proportion of the total number of points wagered on each trial, the proportion of trials on which the more likely outcome was selected (quality of decision making), an average of the proportion wagered on trials when the participant selected the more likely result (rational bets), risk adjustment assessed by variation in the amount wagered in response to the ratio of red to blue boxes, and an index of delay aversion reflected in making higher bets when the amount to bet is presented in descending order rather than in ascending order.
- 'cgt_delay_aversion'  Delay aversion is a measure of the particicipant's preference for making higher bets when the amount to wager is presented in descending (rather than ascending) order.
- 'cgt_deliberation_time' Deliberation time is the time taken to select an amount to wager.
- 'cgt_quality_of_decision_making' Proportion of trials on which the more likely outcome was selected
- 'cgt_overall_proportion_bet' The average of the proportion of the participant's current kitty that is wagered on each trial
- 'cgt_risk_adjustment' The degree to which the proportion bet is affected by the ratio of red:blue boxes on every trial, accounting for the overall proportion bet. (N.B. This measure is not available at baseline for 1 site)
- 'cgt_risk_taking' The average of the proportion of the participant's current kitty that is wagered on each trial, but only for those trials on which they selected the more likely option.