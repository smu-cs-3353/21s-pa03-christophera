# Project 03 - Community Detection

## Overview

Community detection in networks is a popular area of graph theory.  The Girvan-Newman algorithm is a relatively new option for performing community detection in networks/graphs.  

This project implements the Girvan-Newman algorithm using the Boost Graph Library tp detect sub communities within a given graph.  You can find more about Girvan-Newman in the original paper and also some other resources:
* Girvan M. and Newman M. E. J., [Community structure in social and biological networks](https://www.pnas.org/content/99/12/7821), Proc. Natl. Acad. Sci. USA 99, 7821–7826 (2002)
* Chapter 10 Section 2 of [Mining Massive Data Sets](http://infolab.stanford.edu/~ullman/mmds/book0n.pdf)
* [Social and Information Network Analysis Course Slide Deck 14](http://snap.stanford.edu/class/cs224w-2010/slides/14-communities_annot.pdf)

## Requirements 

This project requires two arguments to run:
  1. input file name
  2. output file Name

## Output
Outputs the communities detected into a separate file. 
Each community should have a heading `Community X` followed by a list of nodes that are part of that community each on a separate line. 

You can use the [NetworkX](https://networkx.org/), a Python network analysis library to generate graphs to use for testing.  Take a look at the Graph Generator functionality here The library also has 4 social network graphs already included.  You'll just need to write them out to a file in the format below so your program can read them. 

Example network input:

```text
9
A - B
A - C
B - C
B - D
D - G
D - F
D - E
G - F
E - F
```

