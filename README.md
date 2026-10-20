# Repository Coverage

[Full report](https://htmlpreview.github.io/?https://github.com/scholl-lab/variantcentrifuge/blob/python-coverage-comment-action-data/htmlcov/index.html)

| Name                                                          |    Stmts |     Miss |      Cover |   Missing |
|-------------------------------------------------------------- | -------: | -------: | ---------: | --------: |
| variantcentrifuge/\_\_init\_\_.py                             |        1 |        0 |    100.00% |           |
| variantcentrifuge/\_davies\_build.py                          |        6 |        6 |      0.00% |     11-34 |
| variantcentrifuge/analysis\_profiles.py                       |      138 |        7 |     94.93% |137-139, 263-266 |
| variantcentrifuge/analyze\_variants.py                        |      147 |       71 |     51.70% |71-72, 95-102, 105-170, 206-207, 211-216, 221-227, 230-233, 262, 279-283 |
| variantcentrifuge/annotator.py                                |      174 |       14 |     91.95% |31-34, 181-183, 285, 406, 413, 422, 428, 438, 442 |
| variantcentrifuge/association/\_\_init\_\_.py                 |        4 |        0 |    100.00% |           |
| variantcentrifuge/association/backends/\_\_init\_\_.py        |       14 |        0 |    100.00% |           |
| variantcentrifuge/association/backends/base.py                |       15 |        0 |    100.00% |           |
| variantcentrifuge/association/backends/coast\_python.py       |      162 |       16 |     90.12% |173, 177, 184-186, 201-203, 211-213, 338, 393, 507-509 |
| variantcentrifuge/association/backends/davies.py              |      142 |       13 |     90.85% |86-91, 129-131, 262-263, 307, 429-443 |
| variantcentrifuge/association/backends/python\_backend.py     |      142 |       17 |     88.03% |97-98, 105-106, 113-119, 124, 283, 292, 371-372, 510-511 |
| variantcentrifuge/association/backends/r\_backend.py          |      183 |       35 |     80.87% |53-72, 226-232, 238-239, 267, 273, 437-438, 441, 482-483, 490-491, 532-533, 564 |
| variantcentrifuge/association/backends/skat\_math.py          |      228 |       20 |     91.23% |70, 101, 176, 184-186, 238, 263, 270, 300, 312, 383, 386, 458, 496-501 |
| variantcentrifuge/association/base.py                         |       84 |        0 |    100.00% |           |
| variantcentrifuge/association/correction.py                   |       83 |        9 |     89.16% |26-27, 62-63, 224-225, 238, 252-253 |
| variantcentrifuge/association/covariates.py                   |       48 |        6 |     87.50% |     81-86 |
| variantcentrifuge/association/data/\_\_init\_\_.py            |        0 |        0 |    100.00% |           |
| variantcentrifuge/association/diagnostics.py                  |      173 |       21 |     87.86% |275-276, 281-317, 406, 437 |
| variantcentrifuge/association/engine.py                       |      206 |        6 |     97.09% |200-201, 241, 267, 561-562 |
| variantcentrifuge/association/engine\_workers.py              |       70 |        9 |     87.14% |49-55, 91, 104 |
| variantcentrifuge/association/genotype\_matrix.py             |      162 |       14 |     91.36% |126-133, 211-213, 216, 222-223 |
| variantcentrifuge/association/meta/\_\_init\_\_.py            |        4 |        0 |    100.00% |           |
| variantcentrifuge/association/meta/align.py                   |       50 |        3 |     94.00% |52, 101-102 |
| variantcentrifuge/association/meta/artifact.py                |      144 |       11 |     92.36% |108-110, 143, 146, 148, 171, 320, 326, 329, 355 |
| variantcentrifuge/association/meta/cli.py                     |       28 |        3 |     89.29% |     77-79 |
| variantcentrifuge/association/meta/combine.py                 |       64 |        3 |     95.31% |63, 80, 136 |
| variantcentrifuge/association/meta/heterogeneity.py           |       45 |        3 |     93.33% | 39-40, 82 |
| variantcentrifuge/association/meta/manifest.py                |       62 |        6 |     90.32% |54, 67, 81-82, 97, 114 |
| variantcentrifuge/association/meta/runner.py                  |       63 |        3 |     95.24% |75, 79, 115 |
| variantcentrifuge/association/pca.py                          |      116 |        7 |     93.97% |71, 143, 151-152, 217, 223-224 |
| variantcentrifuge/association/weights.py                      |      198 |        6 |     96.97% |97, 106-107, 151, 161, 404 |
| variantcentrifuge/build\_pm5\_lookup.py                       |      164 |        3 |     98.17% |74, 114, 407 |
| variantcentrifuge/checkpoint.py                               |       84 |        8 |     90.48% |71-74, 94, 102-104 |
| variantcentrifuge/checkpoint\_helpers.py                      |       58 |        0 |    100.00% |           |
| variantcentrifuge/checkpoint\_models.py                       |      106 |        7 |     93.40% |29, 62-65, 82, 110 |
| variantcentrifuge/checkpoint\_state.py                        |      241 |       70 |     70.95% |80-84, 105-107, 127-132, 142, 205-223, 255-261, 288-291, 297-302, 326, 334-338, 342, 346-350, 354-360, 364-383, 409-410, 419, 422, 427, 431-436 |
| variantcentrifuge/cli.py                                      |      205 |       10 |     95.12% |286-288, 297, 300-302, 315, 451, 453 |
| variantcentrifuge/cli\_analysis\_groups.py                    |       98 |        0 |    100.00% |           |
| variantcentrifuge/cli\_config.py                              |      329 |       19 |     94.22% |42, 45-46, 160-161, 268, 275, 295, 304, 313-318, 332, 460, 466, 497 |
| variantcentrifuge/cli\_fields.py                              |       13 |        0 |    100.00% |           |
| variantcentrifuge/cli\_memory.py                              |       14 |        3 |     78.57% |     25-27 |
| variantcentrifuge/cli\_parser\_groups.py                      |       65 |        0 |    100.00% |           |
| variantcentrifuge/cli\_subcommands.py                         |      257 |       23 |     91.05% |139-141, 250, 254, 328, 348, 429-431, 440, 459-471 |
| variantcentrifuge/cli\_sync.py                                |      165 |        0 |    100.00% |           |
| variantcentrifuge/clinvar\_pm5.py                             |       79 |        0 |    100.00% |           |
| variantcentrifuge/config.py                                   |       52 |        2 |     96.15% |    51, 60 |
| variantcentrifuge/converter.py                                |      160 |        8 |     95.00% |120, 248-250, 257, 278-280 |
| variantcentrifuge/converter\_helpers.py                       |       76 |        9 |     88.16% |32, 46-48, 52, 77-79, 126 |
| variantcentrifuge/dataframe\_optimizer.py                     |      141 |        5 |     96.45% |151, 342-344, 349 |
| variantcentrifuge/display\_utils.py                           |      125 |       81 |     35.20% |29-90, 164-197, 209-239 |
| variantcentrifuge/extractor.py                                |      202 |       22 |     89.11% |149-151, 248, 255, 292, 376-377, 394-395, 410-420, 428, 454-455, 561-564 |
| variantcentrifuge/feature\_loaders.py                         |      133 |       16 |     87.97% |93-94, 116-117, 128-129, 173-174, 204-205, 215-216, 230, 234, 253-254 |
| variantcentrifuge/field\_profile.py                           |       48 |        0 |    100.00% |           |
| variantcentrifuge/filters.py                                  |      246 |       24 |     90.24% |69-70, 77-78, 126-127, 133-134, 227-229, 240-241, 260, 263, 275, 287, 397-398, 502-504, 508-509 |
| variantcentrifuge/from\_tsv.py                                |      197 |       16 |     91.88% |90-94, 123, 127, 198, 201-204, 241-242, 245, 281 |
| variantcentrifuge/gene\_bed.py                                |      139 |       36 |     74.10% |52, 55, 58, 66-67, 122-154, 208-209, 274, 288, 290 |
| variantcentrifuge/gene\_burden.py                             |      263 |       27 |     89.73% |16-17, 24-26, 30-31, 82-84, 175, 178, 183-184, 211, 314-315, 356-361, 374-379, 414-419, 447 |
| variantcentrifuge/generate\_html\_report.py                   |      145 |        4 |     97.24% |23, 28, 40, 244 |
| variantcentrifuge/generate\_igv\_report.py                    |      155 |       13 |     91.61% |106-108, 130, 133, 242-243, 258, 266, 289, 302, 310-311 |
| variantcentrifuge/genotype\_utils.py                          |       63 |        6 |     90.48% |38, 44-45, 78, 93, 115 |
| variantcentrifuge/gt\_reconstruction.py                       |       99 |        5 |     94.95% |74-75, 192, 211-215 |
| variantcentrifuge/helpers.py                                  |      216 |      102 |     52.78% |77-89, 119, 123-124, 167, 235-236, 258-259, 265-272, 324-330, 347-348, 376, 379-381, 398-425, 459-495, 513-520, 539-564 |
| variantcentrifuge/inheritance/\_\_init\_\_.py                 |        2 |        0 |    100.00% |           |
| variantcentrifuge/inheritance/analyzer.py                     |      209 |       12 |     94.26% |106-107, 132, 175, 186, 369, 392-394, 489-490, 513-514 |
| variantcentrifuge/inheritance/comp\_het\_vectorized.py        |      201 |       12 |     94.03% |157, 163, 177, 273, 397, 496-499, 563, 567, 572 |
| variantcentrifuge/inheritance/deducer.py                      |      246 |       53 |     78.46% |49, 54, 167, 186, 214, 218, 280, 300, 315, 325, 346, 363-369, 394, 400, 405, 443, 468, 493, 522-553, 581-597 |
| variantcentrifuge/inheritance/non\_autosomal\_deducer.py      |      134 |       16 |     88.06% |183, 217-218, 224-225, 231, 275-279, 326-332 |
| variantcentrifuge/inheritance/parallel\_analyzer.py           |      157 |       17 |     89.17% |45-54, 101, 209-210, 307, 318, 337-339, 351, 362 |
| variantcentrifuge/inheritance/prioritizer.py                  |       44 |        3 |     93.18% |80, 120, 132 |
| variantcentrifuge/inheritance/segregation\_checker.py         |      172 |       20 |     88.37% |70, 89, 92, 118, 122, 130, 140, 175, 179, 202, 228, 238, 263, 275, 278-282, 304, 336, 339 |
| variantcentrifuge/inheritance/vectorized\_deducer.py          |      185 |        7 |     96.22% |94, 118, 122, 311, 340, 367, 527 |
| variantcentrifuge/interactive\_resume.py                      |      189 |      177 |      6.35% |35-134, 159-184, 189-210, 217-266, 273-324, 346-351 |
| variantcentrifuge/links.py                                    |       75 |        7 |     90.67% |134, 166, 169, 187-188, 206-207 |
| variantcentrifuge/memory/\_\_init\_\_.py                      |        2 |        0 |    100.00% |           |
| variantcentrifuge/memory/resource\_manager.py                 |      155 |       17 |     89.03% |129-130, 166-167, 212-225, 237-239 |
| variantcentrifuge/ped\_reader.py                              |       43 |        3 |     93.02% | 43, 63-64 |
| variantcentrifuge/phenotype.py                                |      133 |       20 |     84.96% |69, 77, 79, 87, 123-128, 155-162, 258-261, 308-310 |
| variantcentrifuge/phenotype\_helpers.py                       |      209 |       24 |     88.52% |59-60, 68-69, 80-83, 102-105, 189, 191, 210, 217, 225, 236-238, 251-253, 325 |
| variantcentrifuge/pipeline.py                                 |      208 |       32 |     84.62% |133, 189, 198, 258, 281, 294, 355, 402-406, 421-427, 491, 509-531, 554-557, 572-573 |
| variantcentrifuge/pipeline\_core/\_\_init\_\_.py              |        6 |        0 |    100.00% |           |
| variantcentrifuge/pipeline\_core/context.py                   |      149 |       22 |     85.23% |46, 55-57, 61, 290-292, 308, 310, 312, 314, 316, 318, 320, 324, 326, 332, 334, 350, 354, 364 |
| variantcentrifuge/pipeline\_core/error\_handling.py           |      176 |       23 |     86.93% |90, 192, 197, 223, 240, 253, 272-277, 288, 302-306, 350-351, 442-443, 478-479, 488-489 |
| variantcentrifuge/pipeline\_core/runner.py                    |      191 |       20 |     89.53% |17-18, 164-166, 230, 242, 246, 250, 276, 422-428, 474-483 |
| variantcentrifuge/pipeline\_core/runner\_planner.py           |      196 |       52 |     73.47% |87-98, 103-124, 143, 151, 175-180, 191, 194, 200, 209, 212, 218-235, 242-252 |
| variantcentrifuge/pipeline\_core/stage.py                     |      102 |       12 |     88.24% |173-180, 197-198, 232-235, 246-247 |
| variantcentrifuge/pipeline\_core/workspace.py                 |       62 |        3 |     95.16% |187-188, 221 |
| variantcentrifuge/pipeline\_validation.py                     |       32 |        0 |    100.00% |           |
| variantcentrifuge/pseudonymizer.py                            |      214 |       16 |     92.52% |235, 251, 257, 262, 267-273, 278, 474-475, 479-482 |
| variantcentrifuge/run\_summary.py                             |      161 |       12 |     92.55% |33, 47, 49, 54, 85-86, 114-115, 211-214 |
| variantcentrifuge/sample\_sheet.py                            |      157 |        5 |     96.82% |90, 247-250 |
| variantcentrifuge/scoring.py                                  |      267 |       13 |     95.13% |36, 38, 43, 46, 60-61, 88, 244, 276, 308-309, 416, 421 |
| variantcentrifuge/scoring\_requirements.py                    |       26 |        0 |    100.00% |           |
| variantcentrifuge/stages/\_\_init\_\_.py                      |        5 |        0 |    100.00% |           |
| variantcentrifuge/stages/analysis\_stage\_helpers.py          |       84 |        4 |     95.24% |153, 175-177 |
| variantcentrifuge/stages/analysis\_stages.py                  |       22 |        0 |    100.00% |           |
| variantcentrifuge/stages/association\_analysis\_stages.py     |      284 |       77 |     72.89% |100-122, 162-163, 165-166, 186-189, 214-215, 234-272, 296, 319-323, 329-330, 337, 345-353, 360-375, 404-410, 420-421, 424-425, 439-441, 447-448, 512-517, 538, 542-546 |
| variantcentrifuge/stages/bed\_stages.py                       |      133 |       21 |     84.21% |86, 93-109, 144-159, 166-169, 254-256 |
| variantcentrifuge/stages/chunk\_helpers.py                    |      257 |       65 |     74.71% |31, 36-37, 42-45, 56-58, 79, 81, 92-97, 109, 141-164, 180-184, 221-242, 304-306, 383-392, 414, 431, 441, 452 |
| variantcentrifuge/stages/chunked\_analysis\_stages.py         |      261 |       54 |     79.31% |50, 79-81, 114-128, 156-175, 198, 200, 202, 204, 208, 222-223, 274-279, 385, 395-397, 435-436, 439, 443-445, 454-455, 479, 487-499, 502 |
| variantcentrifuge/stages/data\_loading\_stages.py             |      174 |       53 |     69.54% |26, 56, 90, 94, 106-121, 162, 176-177, 179-180, 186-190, 193-200, 205-207, 213-225, 229, 253, 257, 266-267, 271-272, 293-304, 308 |
| variantcentrifuge/stages/gene\_burden\_stages.py              |      139 |       47 |     66.19% |44-66, 75-85, 98, 133-134, 160, 170-172, 176-179, 182-183, 248, 252-256 |
| variantcentrifuge/stages/genotype\_matrix\_builder.py         |       83 |       11 |     86.75% |61-62, 64, 66, 68, 110-113, 137, 141 |
| variantcentrifuge/stages/inheritance\_and\_scoring\_stages.py |      146 |       19 |     86.99% |49-50, 55-70, 206-207, 213-216, 220, 261-262, 266-267 |
| variantcentrifuge/stages/metadata\_helpers.py                 |       84 |        4 |     95.24% |42, 72, 126, 128 |
| variantcentrifuge/stages/output\_stages.py                    |      292 |       45 |     84.59% |71, 81-82, 86-87, 103, 107-108, 182, 217-218, 221-223, 255-256, 274-275, 286-287, 292-293, 302-303, 308-309, 318-319, 365-366, 383-384, 400, 447-448, 452-453, 495-496, 528-529, 563-566 |
| variantcentrifuge/stages/parallel\_processing\_stage.py       |      340 |       26 |     92.35% |113, 123-126, 172, 204-208, 212, 228, 240, 354, 384, 479-481, 532-536, 569, 593 |
| variantcentrifuge/stages/phenotype\_case\_control\_stages.py  |      285 |       42 |     85.26% |23, 49-52, 54-57, 101-102, 115-117, 128-131, 157, 159, 198, 200, 272, 300-305, 314-315, 350-353, 361-362, 382-384, 395-428, 471-473 |
| variantcentrifuge/stages/post\_processing\_stages.py          |      307 |       64 |     79.15% |49, 53-54, 69-71, 75, 89, 118-120, 126, 130-147, 156, 181, 197-198, 210-214, 218-220, 239, 302-303, 319-320, 324-329, 332-333, 390, 395-396, 450, 455, 464, 477-482, 486-488, 573 |
| variantcentrifuge/stages/processing\_stages.py                |      257 |       29 |     88.72% |95, 100, 130-136, 140, 149-151, 252-253, 259-260, 310-311, 336-337, 341-342, 411, 486-487, 490, 541-543 |
| variantcentrifuge/stages/sample\_and\_pedigree\_stages.py     |      209 |       59 |     71.77% |52, 84-85, 87-91, 113-120, 132, 148, 184, 263-271, 305-308, 314-371 |
| variantcentrifuge/stages/sample\_column\_helpers.py           |      151 |       18 |     88.08% |43-44, 48, 51-52, 161-165, 201-202, 231-236, 339, 347-351, 389 |
| variantcentrifuge/stages/setup\_stages.py                     |      294 |       74 |     74.83% |49-64, 106, 130-131, 142, 174-182, 185-187, 191-199, 203-205, 214-218, 221-227, 254, 288, 297-303, 333-335, 349, 365, 391-392, 413-414, 421-423, 437, 461, 473, 481, 490, 514, 522, 530-535 |
| variantcentrifuge/stages/stage\_registry.py                   |      183 |       70 |     61.75% |32, 34, 36, 70, 77-80, 84, 92-93, 95, 141, 156-157, 172-176, 186, 201-202, 212, 227, 242-256, 275-310, 329-341, 351-365 |
| variantcentrifuge/stages/stage\_utils.py                      |       20 |       20 |      0.00% |     9-101 |
| variantcentrifuge/stages/statistics\_and\_clinvar\_stages.py  |      122 |       47 |     61.48% |47-60, 71-72, 76, 106-108, 111-113, 117-119, 123-124, 213-241 |
| variantcentrifuge/stages/tsv\_output\_stages.py               |      252 |       30 |     88.10% |76-80, 106, 122-123, 148-155, 174, 193, 260-264, 270-273, 392-393, 417-418, 431, 447 |
| variantcentrifuge/stages/variant\_analysis\_stages.py         |       93 |        6 |     93.55% |124, 129, 165-168 |
| variantcentrifuge/stats.py                                    |       87 |        5 |     94.25% |115, 149-150, 176-177 |
| variantcentrifuge/stats\_engine.py                            |      199 |       42 |     78.89% |228, 248, 253, 358-361, 388-389, 426-428, 432-435, 444-467, 504-506, 519-522, 545-549 |
| variantcentrifuge/transcript\_filter.py                       |      222 |       14 |     93.69% |29, 45, 67, 112, 136, 155, 167, 184, 187, 209, 215, 222, 250, 308 |
| variantcentrifuge/tsv\_filters.py                             |      207 |      162 |     21.74% |95-304, 345-348, 353-370, 380-386, 393-395, 415 |
| variantcentrifuge/utils.py                                    |      212 |       51 |     75.94% |39-42, 66-68, 74, 119-120, 136-139, 223-226, 239-242, 270-271, 382, 388, 464, 467, 474-476, 483-487, 541-547, 570-584 |
| variantcentrifuge/validators.py                               |      300 |       29 |     90.33% |69, 93, 95, 123, 125, 133, 137-142, 152, 154, 163, 168, 170, 212-214, 255, 263, 265, 273, 279, 285, 291, 298, 302 |
| variantcentrifuge/vcf\_eff\_one\_per\_line.py                 |      283 |       71 |     74.91% |104, 148, 173, 178-179, 196-210, 228-234, 251, 255, 259, 266-272, 293, 302, 345, 393-447, 452-472 |
| variantcentrifuge/vcf\_header\_parser.py                      |      171 |       18 |     89.47% |153, 164, 199, 204, 231-245, 380, 392-393, 410 |
| variantcentrifuge/version.py                                  |        1 |        0 |    100.00% |           |
| **TOTAL**                                                     | **17556** | **2703** | **84.60%** |           |


## Setup coverage badge

Below are examples of the badges you can use in your main branch `README` file.

### Direct image

[![Coverage badge](https://raw.githubusercontent.com/scholl-lab/variantcentrifuge/python-coverage-comment-action-data/badge.svg)](https://htmlpreview.github.io/?https://github.com/scholl-lab/variantcentrifuge/blob/python-coverage-comment-action-data/htmlcov/index.html)

This is the one to use if your repository is private or if you don't want to customize anything.

### [Shields.io](https://shields.io) Json Endpoint

[![Coverage badge](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/scholl-lab/variantcentrifuge/python-coverage-comment-action-data/endpoint.json)](https://htmlpreview.github.io/?https://github.com/scholl-lab/variantcentrifuge/blob/python-coverage-comment-action-data/htmlcov/index.html)

Using this one will allow you to [customize](https://shields.io/endpoint) the look of your badge.
It won't work with private repositories. It won't be refreshed more than once per five minutes.

### [Shields.io](https://shields.io) Dynamic Badge

[![Coverage badge](https://img.shields.io/badge/dynamic/json?color=brightgreen&label=coverage&query=%24.message&url=https%3A%2F%2Fraw.githubusercontent.com%2Fscholl-lab%2Fvariantcentrifuge%2Fpython-coverage-comment-action-data%2Fendpoint.json)](https://htmlpreview.github.io/?https://github.com/scholl-lab/variantcentrifuge/blob/python-coverage-comment-action-data/htmlcov/index.html)

This one will always be the same color. It won't work for private repos. I'm not even sure why we included it.

## What is that?

This branch is part of the
[python-coverage-comment-action](https://github.com/marketplace/actions/python-coverage-comment)
GitHub Action. All the files in this branch are automatically generated and may be
overwritten at any moment.