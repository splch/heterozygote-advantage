# Heterozygote Advantage Simulation
Simulations of populations with varying levels of natural selection to observe responses to heterozygous percentages.

## Hypothesis
The rate of natural selection is directly proportional to the frequency of heterozygous genes.

## Installation

[![Run on Repl.it](https://replit.com/badge/github/splch/heterozygote-advantage)](https://replit.com/github/splch/heterozygote-advantage)

```bash
pip install numpy scipy matplotlib pysimplegui
```

## Usage
<img src="Images/GUI.png" height="400px">

To run application:
```bash
python hetadv.py
```

Then simply set:
  number of genes,
  birth rate,
  generations,
  and iterations

*Keep in mind this program slows exponentially for more generations (especially above 20).*

## Results
<img src="Images/Figure_1.png" height="400px">

Here we can see a strong correlation between increasing natural selection and frequency heterozygous genes.
