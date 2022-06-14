# HUST_RTS_SOLVER

本求解器为PACE2022的参赛之作。（[PACE2022官网](https://pacechallenge.org/2022/)）
The solver is a submission for PACE2022. ([PACE2022 Official Website](https://pacechallenge.org/2022/))

## 算法说明 Algorithm Description

6月15日前补充pdf链接。
Add pdf link by June 15.

## 编译 Compilation

```shell
mkdir build
cd build
cmake ..
make
```

编译生成符合[optil](https://www.optil.io/optilion/problem/3198)要求的程序。
Generate program that meet the requirements of [optil](https://www.optil.io/optilion/problem/3198).

## 执行 Execution

`.\PACE`
接受命令行输入，将结果输出至命令行中。
Take input from command lines and output to command lines.

## 输入格式 Input Format

第一行第一个整数为点的个数N
接下来N行中的第i行表示从点i起始的边的终点编号，以空格分开
点的编号为1~N
例：

> 5
> 2 3 4
> 3 5
>
> 1 2
> 4

## 输出格式 Output Format

一系列整数，表示找到的最小反馈集所包含的点编号

## 参考文献 Reference

- Philippe Galinier, Eunice Lemamou, and Mohamed Wassim Bouzidi. Applying local search to the feedback vertex set problem. Journal of Heuristics, 19(5):797–818, oct 2013. doi: 10.1007/s10732-013-9224-z.
- Hen-Ming Lin and Jing-Yang Jou. On computing the minimum feedback vertex set of a directed graph by contraction operations. IEEE Transactions on Computer-Aided Design of Integrated Circuits and Systems, 19(3):295–307, mar 2000. doi:10.1109/43.833199.
- Levy, H., Low, D.W.: A contraction algorithm for finding small cycle cutsets. J. Algorithms 9(4), 470–493 (1988).