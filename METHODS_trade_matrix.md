# Generating the transition matrix
The transition matrix is necessary for the simulation to run. It defines the probability of a host to move between discrete locations, where each row and column sums to 1. The question centers on sampling bias and uneven trade flows, and this should be reflected in the matrix. The current approach is to:

1. Convert trade quantities from metric tonnes to estimated numbers of pigs using assumed average weights for each swine category.
2. Identify Denmark’s five largest export destinations over a three-year trade period.
3. Restrict the trade network to Denmark and these five partner countries.
4. Aggregate trade between each country pair across the three years and convert three-year totals to average daily trade.
5. Obtain herd sizes using 2024 November–December data for Live swine, domestic species.
6. Calculate per-pig movement rates as daily exports divided by the exporter’s herd size.
7. Normalize movement rates within each exporter to obtain destination probabilities conditional on a movement occurring, by dividing the route's rate by the sum of the exporters rate.
   This answers the question: if a movement from country X occurs, what is the probability that this movement is to country Y?
8. Populate matrix

# Assumptions
Conversion from tonnes to animals
    The assumed weights.
    Why those weights were selected.

Temporal conversion
    Three-year totals → average daily movement.
    Why a daily timestep is being used.

# Data sources
Trade data
    Source.
    Years included.
    What quantity represents.
    Which commodity categories are included.
Conversion from tonnes to animals
    The assumed weights.
    Why those weights were selected.
    The limitations of using fixed mean weights.
Temporal conversion
    Three-year totals → average daily movement.
    Why a daily timestep is being used.
Herd population
    Source.
    Why 2024 is used.
    Why November–December is used.
    Why Live swine, domestic species is used rather than summing the four categories.
Movement rate
Structure matrix

Important distinction

    Explicitly document:

    The structure matrix describes the destination distribution of movements and does not encode the absolute frequency of movement. Absolute movement rates are retained separately.

Limitations

    This is particularly valuable for your eventual paper/thesis. For example:

    fixed animal weights;
    annual/three-year trade data converted to daily averages;
    2024 herd size used as a proxy for the three-year period;
    selected-country network rather than complete global trade;
    excluded trade routes outside the simulation network;
    trade data represent reported movements rather than necessarily actual biological movements;
    potential uncertainty in the conversion from tonnes to numbers of animals.
