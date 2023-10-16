# Scalable Approximate Butterfly and Bi-triangle Counting for Large Bipartite Networks

This repository implements one-sided weighted local sampling method that first proposed in our paper, as well as other existing approximate sampling methods for comparison. As a summary, we implement the following method:

-  approximate butterfly counting
   -  `fast edge sampling`(The original version is from https://github.com/beginner1010/butterfly-counting, but has been refactored into this project)
   -  *`one-sided pair sampling`
   -  *`weighted one-sided pair sampling`
-  approximate bi-triangle counting
   -  `vertex sampling`
   -  `edge sampling`
   -  `wedge sampling`
   -  *`one-sided triple sampling`
   -  *`weighted one-sided triple sampling`

*methods proposed in our paper

## Compile

```shell
make
```

## File Description

-  `app.cpp`: Program entrance and argument parsing.
-  `bi_graph.*`: Reading and processing of bipartite graphs, including local sampling methods.
-  `exact_alg.*`: Implementation of exact algorithms for butterfly counting and bi-triangle counting.
-  `tracker.*`: Tracking sampling algorithms, recording the time required for each number of samples and the corresponding error rate.
-  `gt_manager.*`: Managing ground truth.
-  `util.*`: Utility functions.
-  `def.h`: Definition of constants and type aliases, etc.

## Parameters

-  `-alias [STRING]`: Dataset alias
   -  The program will read the graph from `DATASET_ROOT/out.<alias>`
   -  `DATASET_ROOT` is defined in `def.h`, please motify it before compiling
-  `-me [bfc|btc]`: Patterns that require quantity counting
   -  `-me bfc` for butterfly counting
   -  `-me btc` for bi-triangle counting
-  `-alg [1~5]`: Number of the algorithm used
   -  `-alg 1` for vertex sampling
   -  `-alg 2` for (fast) edge sampling
   -  `-alg 3` for wedge sampling
   -  `-alg 4` for pair/triple sampling for butterfly counting/bi-triangle counting respectively  
   -  `-alg 5` for weighted pair/triple sampling for butterfly counting/bi-triangle counting respectively  
-  `-rnd [INT]`: Number of running rounds (each round runs independently)
-  `-side [l|r|b]`: Sampling mode, **not valid** for `(fast) edge sampling`
   -  `l` means sampling from the left side
   -  `r` means sampling from the right side
   -  `b` means sampling from both sides

## Input & Output

The graph must be stored in `DATASET_ROOT/out.<alias>` in the **edge list** format and `DATASET_ROOT` is defined in `def.h`. The vertices on both sides are numbered independently. For example:

```
1 1
1 2
1 3
2 3
```

The result will be outputed to `RESULT_ROOT/<me>_<alg_name>.<alias><rnd>[_<side>]`, and `RESULT_ROOT` is defined in `def.h`. The output is four columns of data: time, number of samples, error, and estimate. To get error, the ground truth needs to be prepared in advance in `GT_ROOT/gt.<alias>`, and `GT_ROOT` is defined in `def.h`. The file `gt.<alias>` is three columns of data: name, value and time (optional). 

For example, if we want to get error directly when testing the approximation algorithm for butterfly counting on the dataset `orkut`, we need to add a line `bfc 22131701213295 340.95` to the file `GT_ROOT/gt.orkut`. 

