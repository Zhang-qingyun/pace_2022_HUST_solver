# HUST_RTS_SOLVER

The solver is a submission for PACE2022. ([PACE2022 Official Website](https://pacechallenge.org/2022/))

## Algorithm Description

See `Solver_description-HUST_RTS_SOLVER.pdf`.

[![10.5281/zenodo.6643002](https://zenodo.org/badge/503428440.svg)](https://doi.org/10.5281/zenodo.6643002)

```

```

## Compilation

```shell
mkdir build
cd build
cmake ..
make
Generate program that meet the requirements of [optil](https://www.optil.io/optilion/problem/3198).
```

## Execution

`.\PACE22_FVSP_HUST_NOSCP`
Take input from command lines and output to command lines.

## Input Format

> graphs are considered to be directed and simple, i.e. parallel edges and self-loops are forbidden. For example, a graph with four vertices and five directed edges can be defined as follows:
>
> ```
> 4 5 0 
> 2 3 
> 3
> 4 
> 1
> ```

## Output Format

Every line must be in the form `u` followed by the new line character `\n`, where `u` represents a node contained in your feedback vertex set. Here is an example:

```
 1 
 2
```

## Reference

- Philippe Galinier, Eunice Lemamou, and Mohamed Wassim Bouzidi. Applying local search to the feedback vertex set problem. Journal of Heuristics, 19(5):797–818, oct 2013. doi: 10.1007/s10732-013-9224-z.
- Hen-Ming Lin and Jing-Yang Jou. On computing the minimum feedback vertex set of a directed graph by contraction operations. IEEE Transactions on Computer-Aided Design of Integrated Circuits and Systems, 19(3):295–307, mar 2000. doi:10.1109/43.833199.
- Levy, H., Low, D.W.: A contraction algorithm for finding small cycle cutsets. J. Algorithms 9(4), 470–493 (1988).