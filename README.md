## Efficiently merging two r-indexes

#### Abstract
Large sequencing projects, such as GenomeTrakr and MetaSub, are updated frequently(sometimes daily, in the case of GenomeTrakr) with new data.  Therefore, it is imperativethat any data structure indexing such data supports efficient updates.  Toward this goal, Bannai  et  al.   (TCS,  2020)  proposed  a  data  structure  named dynamic r-index  which  is suitable for large genome collections and supports incremental construction; however, it is still not powerful enough to support substantial updates. Here, we develop a novel algorithm for updating ther-index, which we refer to as `rimerge`.  Fundamental to our algorithm is the combination of the basics of the dynamic r-index with a known algorithm for merging Burrows-Wheeler Transforms (BWTs).  As a result, rimerge is capable of performing batch updates in a manner that exploits parallelism while keeping the memory overhead small. We compare our method to the dynamic r-index of Bannai et al.  using two different datasets, and show that `rimerge` is between 1.88 to 5.34 times faster on reasonably large inputs.


#### Donwload and run

In order to compile and run `rimerge` run the following:

```bash
git clone https://github.com/marco-oliva/rimerge
cd rimerge
mkdir build && cd build
cmake ..
make 
make install
```
This will create a `bin` directory where you will find `rimerge`


