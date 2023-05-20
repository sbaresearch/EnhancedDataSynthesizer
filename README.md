# Enhanced DataSynthesizer

This is an extended and enhanced version of the [DataSynthesizer](https://github.com/DataResponsibly/DataSynthesizer), which generates synthetic data that simulates a given
dataset. For details on that project and usage, please refer to the official website!

## Extensions

This version extends the original DataSynthesizer by two features:

### Fast building of Bayesian networks
The original DataSynthesizer relies on a greedy approach to find a fitting network structure, which becomes (prohibitedly) expensive with larger numbers of samples, features, or
parent nodes. This enhanced version provides a mode using Genetic Algorithms to perform a heuristic, much faster search for a fitting network structure

### Protection of sensitive attributes
By specifically reducing or even removing correlations between sensitive attributes and other attributes, attacks like attribute inference are reduced.

Details on these modifications can also be found in the paper [Efficient Bayesian Network Construction for Increased Privacy on Synthetic Data](http://dx.doi.org/10.1109/BigData55660.2022.10020936) by Markus Hittmeir, Rudolf Mayer, and Andreas Ekelhart.

To cite that paper, please use the following:

```
@inproceedings{hittmeir_efficient_2022,
  title = {Efficient {Bayesian} {Network} {Construction} for {Increased} {Privacy} on {Synthetic} {Data}},
  booktitle = {2022 {IEEE} {International} {Conference} on {Big} {Data} ({Big} {Data})},
  author = {Markus Hittmeir and Rudolf Mayer and Andreas Ekelhart},
  publisher = {IEEE Computer Society},
  address = {Osaka, Japan},
  month = dec,
  year = {2022},
  doi = {10.1109/BigData55660.2022.10020936},
}

```
