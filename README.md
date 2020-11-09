## Efficiently merging two r-indexes

#### Abstract
AbstractLarge sequencing projects, such as GenomeTrakr and MetaSub, are updated frequently(sometimes daily, in the case of GenomeTrakr) with new data.  Therefore, it is imperativethat any data structure indexing such data supports efficient updates.  Toward this goal,Bannai  et  al.   (TCS,  2020)  proposed  a  data  structure  nameddynamic r-index  which  issuitable for large genome collections and supports incremental construction; however, it isstill not powerful enough to support substantial updates. Here, we develop a novel algorithmfor updating ther-index, which we refer to asrimerge.  Fundamental to our algorithm isthe combination of the basics of thedynamic r-index with a known algorithm for mergingBurrows-Wheeler Transforms (BWTs).  As a result,rimergeis capable of performing batchupdates in a manner that exploits parallelism while keeping the memory overhead small. Wecompare our method to the dynamicr-index of Bannai et al.  Using two different datasets,and show thatrimergeis between 1.88 to 5.34 times faster on reasonably large inputs.


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
This will create a bin directory where you will find `rimerge`


