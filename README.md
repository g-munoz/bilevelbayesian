# Exploiting the polyhedral geometry of stochastic linear bilevel programming

Arxiv paper: https://arxiv.org/abs/2211.02268

usage: julia main.py filename timelimit 

- filename is either a julia file with the instance or the string "knapsack". In the latter is uses "ContinuousKnapsack_instance.jl" to generate the instance
- timelimit specifies the time limit of the exact method

Sample usages:

*julia main.jl JuliaInstances/BOLib2/BardFalk1982Ex2.jl 600*

*julia main.jl knapsack 600*
