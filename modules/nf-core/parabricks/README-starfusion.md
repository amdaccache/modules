# Parabricks Starfusion

`Starfusion` performs fusion calling for RNA-Seq samples, utilizing the STAR-Fusion algorithm. This requires input of a genome resource library, in accordance with the original STAR-Fusion tool, and outputs candidate fusion transcripts. 
alignment, sorting, (optional) marking of duplicates, and (optional) base quality score recalibration (BQSR). There is no option to control the number of threads used with this tool - all available threads on the system are used by default.

For additional considerations see the [tool documentation](https://docs.nvidia.com/clara/parabricks/latest/documentation/tooldocs/man_starfusion.html).
