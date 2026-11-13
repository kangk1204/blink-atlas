# Additional unseen-parent check

The primary protocol is in transfer-registration.json. After the first registered fits, we added one diagnostic to directly test the newly added chemical identities together. The results of this additional test had not been computed when this file was committed.

Hold out Sakabe, Uno 2018 and Chi 2020 simultaneously. Remove all parent identities in their curated panels from the remaining studies, including parents those papers cite as controls. The primary test targets are every exact single-transition absorption parent from the additions that had no eligible label in the previous 26-parent release. Report all per-parent scores, mean absolute errors, bias and win counts. Average per parent. This exploratory test has no role in model selection.

For HMbisTFER the authors give a 7.2 midpoint, a separate absorbance series in Fig. S3c and a fluorescence series in Fig. S3d. Their methods allow either channel and the accessible figure caption does not specify which channel produced 7.2. Keep that value in a separately labeled, more inclusive sensitivity analysis, both as a training observation and as an external test case. Publish the effect of this decision on **both** prediction methods.
