# iaintersect
Based on
[BioWardrobe iaintersect](https://github.com/cincinnati-childrens-hospital/biowardrobe/tree/master/src/iaintersect "BioWardrobe")

The software assigns each peak obtained from MACS2 to a gene and region (upstream, promoter, exon, intron, intergenic).

## Usage
`iaintersect [options] --in=pathToFile --a=pathtoFile --out=pathToFile`

|Option          |Type    |Info                                                |
|----------------|--------|----------------------------------------------------|
|--in            |Required| Input filename with MACS2 peak calling results, xls|
|--out           |Required| Base output file name, tsv                         |
|--a             |Required| Annotation file, tsv                               |
|--log           |Optional| Log filename (default ./logfile_def.log)           |
|--promoter      |Optional| Promoter region around TSS, base pairs             |
|--upstream      |Optional| Upstream region before promoter, base pairs        |
|--sam_ignorechr |Optional| The chromosome to be ignored, string               |
 
