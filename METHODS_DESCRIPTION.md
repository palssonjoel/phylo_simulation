# Generating the transition matrix

Relevant code: 
 - filter_trade_data.R
 - create_trans_matrix.R 
 
The transition matrix is necessary for the simulation to run. It defines the probability of a host to move between discrete locations, where each row and column sums to 1. The question centers on sampling bias and uneven trade flows, and this should be reflected in the matrix. Thus, the structure matrix describes the destination distribution of movements and does not encode the absolute frequency of movement. Absolute movement rates are retained separately.

The current approach is to:

1. Convert trade quantities from metric tonnes to estimated numbers of pigs using assumed average weights for each swine category.
2. Identify Denmark’s five largest export destinations over a three-year trade period.
3. Restrict the trade network to Denmark and these five partner countries.
4. Aggregate trade between each country pair across the three years and convert three-year totals to average daily trade.
5. Obtain herd sizes using 2024 November–December data for Live swine, domestic species.
6. Calculate per-pig movement rates as daily exports divided by the exporter’s herd size.
7. Normalize movement rates within each exporter to obtain destination probabilities conditional on a movement occurring, by dividing the route's rate by the sum of the exporters rate.
   This answers the question: if a movement from country X occurs, what is the probability that this movement is to country Y?
8. Populate matrix

## Assumptions
The trade data showed live animal trades in tonnes, and has three categories: breeding sows; non-breeding pigs over 50 kg; non-breeding pigs under 50 kg. Average pig weights for breeding pigs and non-breeding pigs were set as 250 kg/pig (sows) and 120 kg/pig (non-breeding adult pigs), and for pigs under 50 kg a weight of 30 kg was assumed. These weights are based on averages found online.

The trade data covered three years: 2022-2024. The trade volumes of these years were aggregated, and to get a daily average this was divided by 365 * 3. Daily average is used because the simulation timesteps will be in days.

## Data sources
Trade data: 
    Source.
    Years included.
    What quantity represents.
    Which commodity categories are included.

Herd population
    Source.
    Why 2024 is used.
    Why November–December is used.
    Why Live swine, domestic species is used rather than summing the four categories.

## Limitations
 - fixed animal weights;
 - annual/three-year trade data converted to daily averages;
 - 2024 herd size used as a proxy for the three-year period;
 - selected-country network rather than complete global (or EU) trade. The matrix describes only the probability of hosts moving within these countries, and any movement outside this is effectively ignored. Actual probability would differ, and demands the use of all available trade data (which is possible to do, but deemed unnecessary).
 - excluded trade routes outside the simulation network;
 - potential uncertainty in the conversion from tonnes to numbers of animals.

## Possible improvements
 - Structure code into a function, list of countries can be an argument, allows better scalability. 
 
# Simulating a transmission chain

Relevant code: 
 - simulation.R

## Description 
Using the transition matrix, a transmission chain is simulated using the nosoi package. This is an agent-based, stochastic transmission chain simulator developed by Lquime et al. It uses discrete space which is defined by the countries present in the transition matrix, with one initial infected individual in Denmark. 

Source:
Sebastian Lequime, Paul Bastide, Simon Dellicour, Philippe Lemey & Guy Baele (2020) nosoi: A stochastic agent-based transmission chain simulation framework in R. Methods in Ecology and Evolution 11:1002-1007 doi:10.1111/2041-210X.13422

## Assumptions
 - Constant probability of transmission once a host is infectious
 - Transmission period is same as disease lenght, minus incubation time. In reality, virus' shed beyond this period.
 - All variables are the same for all locations
 
# HKY substituion model
An HKY substition model is applied to the transmission chain produced by nosoi. Transition/transversion rate differences are defined by the ratio kappa. Parts of Layan et al. HKY code has been adapted here.

Logic:
For each host, it checks who they were infected by and pulls their sequence. A substituion rate is then applied to this sequence for the amount of time between when the infector was infected, and when the host in quesion was infected. The logic here is the same as Layan et al. 

## Assumptions
 - Sites evolve independently
 - Time-homogeneity: substitution rate matrix is constant across entire simulation
 - Base frequencies are stationary according to the reference genomes frequencies across the simulation. Drift not allowed.
 - Kappa is fixed and uniform across sites
 - All sites evolve at the same rate.
 - No selection, purely neutral drift.


Source: github.com/mlayan/Sampling_bias












