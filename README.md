# Project 03 - Community Detection

**Due: Monday April 12, 2021 @ 6am**

## Overview

Community detection in networks is a popular area of graph theory.  The Girvan-Newman algorithm is a relatively new option for performing community detection in networks/graphs.  

For this project, you'll implement the Girvan-Newman algorithm using the Boost Graph Library.  You can find more about Girvan-Newman in the original paper and also some other resources:
* Girvan M. and Newman M. E. J., [Community structure in social and biological networks](https://www.pnas.org/content/99/12/7821), Proc. Natl. Acad. Sci. USA 99, 7821–7826 (2002)
* Chapter 10 Section 2 of [Mining Massive Data Sets](http://infolab.stanford.edu/~ullman/mmds/book0n.pdf)
* [Social and Information Network Analysis Course Slide Deck 14](http://snap.stanford.edu/class/cs224w-2010/slides/14-communities_annot.pdf)


What you need to do:

  1. read about Girvan Newman and figure out how the algorithm works. 
  2. create a new project and explore the Boost Graph Library
  3. implement
    1. reading in the network/graph file based on the format below
    2. detect communities
    3. output the communities detected into a separate file.  Each community should have a heading `Community X` followed by a list of nodes that are part of that community each on a separate line. 

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

Important:  You'll need to include a Github actions workflow file so that your the project will build under Actions.