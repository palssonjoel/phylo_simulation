# Generating the transition matrix
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
