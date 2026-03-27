mapping = {
"basic_statistics": "BASIC STATISTICS",
"per_base_sequence_quality": "PER BASE SEQUENCE QUALITY",
"per_tile_sequence_quality": "PER TILE SEQUENCE QUALITY",
"per_sequence_quality_scores": "PER SEQUENCE QUALITY SCORES",
"per_base_sequence_content": "PER BASE SEQUENCE CONTENT",
"per_sequence_gc_content": "PER SEQUENCE GC CONTENT",
"per_base_n_content": "PER BASE N CONTENT",
"sequence_length_distribution": "SEQUENCE LENGTH DISTRIBUTION",
"sequence_duplication_levels": "SEQUENCE DUPLICATION LEVELS",
"overrepresented_sequences": "OVERREPRESENTED SEQUENCES",
"adapter_content": "ADAPTER CONTENT"
}

reasons_mapping = {
"basic_statistics": "BASIC STATISTICS",

"per_base_sequence_quality":
"""A warning is issued if the lower quartile for any base is less than 10, or if the median for any base is less than 25.
This module will raise a failure if the lower quartile for any base is less than 5 or if the median for any base is less than 20.

Common reasons for warnings:
The most common reason for warnings and failures in this module is a general degradation of quality over the duration of long runs. In general sequencing chemistry degrades with increasing read length and for long runs you may find that the general quality of the run falls to a level where a warning or error is triggered.

If the quality of the library falls to a low level then the most common remedy is to perform quality trimming where reads are truncated based on their average quality. For most libraries where this type of degradation has occurred you will often be simultaneously running into the issue of adapter read-through so a combined adapter and quality trimming step is often employed.

Another possibility is that a warn / error is triggered because of a short loss of quality earlier in the run, which then recovers to produce later good quality sequence. This can happen if there is a transient problem with the run (bubbles passing through a flowcell for example). You can normally see this type of error by looking at the per-tile quality plot (if available for your platform). In these cases trimming is not advisable as it will remove later good sequence, but you might want to consider masking bases during subsequent mapping or assembly.

If your library has reads of varying length then you can find a warning or error is triggered from this module because of very low coverage for a given base range. Before committing to any action, check how many sequences were responsible for triggering an error by looking at the sequence length distribution module results.
""",

"per_tile_sequence_quality":
"""This module issues a warning if any tile shows a mean Phred score more than 2 less than the mean for that base across all tiles.
This module will issue a warning if any tile shows a mean Phred score more than 5 less than the mean for that base across all tiles.

Common reasons for warnings:
Whilst warnings in this module can be triggered by individual specific events we have also observed that greater variation in the phred scores attributed to tiles can also appear when a flowcell is generally overloaded. In this case events appear all over the flowcell rather than being confined to a specific area or range of cycles. We would generally ignore errors which mildly affected a small number of tiles for only 1 or 2 cycles, but would pursue larger effects which showed high deviation in scores, or which persisted for several cycles.
""",

"per_sequence_quality_scores": 
"""A warning is raised if the most frequently observed mean quality is below 27 - this equates to a 0.2% error rate.
An error is raised if the most frequently observed mean quality is below 20 - this equates to a 1% error rate.

Common reasons for warnings:
This module is generally fairly robust and errors here usually indicate a general loss of quality within a run. For long runs this may be alleviated through quality trimming. If a bi-modal, or complex distribution is seen then the results should be evaluated in concert with the per-tile qualities (if available) since this might indicate the reason for the loss in quality of a subset of sequences.
""",

"per_base_sequence_content": 
"""Common reasons for warnings:
The most common reason for warnings and failures in this module is a general degradation of quality over the duration of long runs. In general sequencing chemistry degrades with increasing read length and for long runs you may find that the general quality of the run falls to a level where a warning or error is triggered.

If the quality of the library falls to a low level then the most common remedy is to perform quality trimming where reads are truncated based on their average quality. For most libraries where this type of degradation has occurred you will often be simultaneously running into the issue of adapter read-through so a combined adapter and quality trimming step is often employed.

Another possibility is that a warn / error is triggered because of a short loss of quality earlier in the run, which then recovers to produce later good quality sequence. This can happen if there is a transient problem with the run (bubbles passing through a flowcell for example). You can normally see this type of error by looking at the per-tile quality plot (if available for your platform). In these cases trimming is not advisable as it will remove later good sequence, but you might want to consider masking bases during subsequent mapping or assembly.

If your library has reads of varying length then you can find a warning or error is triggered from this module because of very low coverage for a given base range. Before committing to any action, check how many sequences were responsible for triggering an error by looking at the sequence length distribution module results.
""",

"per_sequence_gc_content": 
"""A warning is raised if the sum of the deviations from the normal distribution represents more than 15% of the reads.
This module will indicate a failure if the sum of the deviations from the normal distribution represents more than 30% of the reads.

Common reasons for warnings:
Warnings in this module usually indicate a problem with the library. Sharp peaks on an otherwise smooth distribution are normally the result of a specific contaminant (adapter dimers for example), which may well be picked up by the overrepresented sequences module. Broader peaks may represent contamination with a different species.
""",

"per_base_n_content": 
"""This module raises a warning if any position shows an N content of >5%.
This module will raise an error if any position shows an N content of >20%.

Common reasons for warnings:
The most common reason for the inclusion of significant proportions of Ns is a general loss of quality, so the results of this module should be evaluated in concert with those of the various quality modules. You should check the coverage of a specific bin, since it is possible that the last bin in this analysis could contain very few sequences, and an error could be prematurely triggered in this case.

Another common scenario is the incidence of a high proportions of N at a small number of positions early in the library, against a background of generally good quality. Such deviations can occur when you have very biased sequence composition in the library to the point that base callers can become confused and make poor calls. This type of problem will be apparent when looking at the per-base sequence content results.
""",

"sequence_length_distribution": 
"""This module will raise a warning if all sequences are not the same length.
This module will raise an error if any of the sequences have zero length.

Common reasons for warnings:
For some sequencing platforms it is entirely normal to have different read lengths so warnings here can be ignored.""",

"sequence_duplication_levels": 
"""This module will issue a warning if non-unique sequences make up more than 20% of the total.
This module will issue a error if non-unique sequences make up more than 50% of the total.

Common reasons for warnings:
The underlying assumption of this module is of a diverse unenriched library. Any deviation from this assumption will naturally generate duplicates and can lead to warnings or errors from this module.

In general there are two potential types of duplicate in a library, technical duplicates arising from PCR artefacts, or biological duplicates which are natural collisions where different copies of exactly the same sequence are randomly selected. From a sequence level there is no way to distinguish between these two types and both will be reported as duplicates here.

A warning or error in this module is simply a statement that you have exhausted the diversity in at least part of your library and are re-sequencing the same sequences. In a supposedly diverse library this would suggest that the diversity has been partially or completely exhausted and that you are therefore wasting sequencing capacity. However in some library types you will naturally tend to over-sequence parts of the library and therefore generate duplication and will therefore expect to see warnings or error from this module.

In RNA-Seq libraries sequences from different transcripts will be present at wildly different levels in the starting population. In order to be able to observe lowly expressed transcripts it is therefore common to greatly over-sequence high expressed transcripts, and this will potentially create large set of duplicates. This will result in high overall duplication in this test, and will often produce peaks in the higher duplication bins. This duplication will come from physically connected regions, and an examination of the distribution of duplicates in a specific genomic region will allow the distinction between over-sequencing and general technical duplication, but these distinctions are not possible from raw fastq files. A similar situation can arise in highly enriched ChIP-Seq libraries although the duplication there is less pronounced. Finally, if you have a library where the sequence start points are constrained (a library constructed around restriction sites for example, or an unfragmented small RNA library) then the constrained start sites will generate huge dupliction levels which should not be treated as a problem, nor removed by deduplication. In these types of library you should consider using a system such as random barcoding to allow the distinction of technical and biological duplicates.
""",

"overrepresented_sequences": 
"""This module will issue a warning if any sequence is found to represent more than 0.1% of the total.
This module will issue an error if any sequence is found to represent more than 1% of the total.

Common reasons for warnings:
This module will often be triggered when used to analyse small RNA libraries where sequences are not subjected to random fragmentation, and the same sequence may natrually be present in a significant proportion of the library.
""",

"adapter_content": 
"""This module will issue a warning if any sequence is present in more than 5% of all reads.
This module will issue an error if any sequence is present in more than 10% of all reads.

Common reasons for warnings:
Any library where a reasonable proportion of the insert sizes are shorter than the read length will trigger this module. This does not indicate a problem as such - just that the sequences will need to be adapter trimmed before proceeding with any downstream analysis.
"""
}
