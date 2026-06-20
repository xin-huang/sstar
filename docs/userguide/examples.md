```
sstar2 train --demes examples/HumanNeanderthal_4G21_wo_introgression.yaml --config examples/sstar2.example.config.yaml --output examples/sstar2.example.trained.model.onnx
sstar2 infer --model examples/sstar2.example.trained.model.onnx --config examples/sstar2.example.config.yaml --feat-file examples/sstar2.examples.inference.features.tsv --pred-file examples/sstar2.example.pred.tsv --tract-file examples/sstar2.example.inferred.tracts.bed
sstar2 match --vcf examples/sstar2.example.biallelic.snps.vcf.gz --tgt examples/tgt.samples.list --src examples/nean.samples.list --tract-file examples/sstar2.example.inferred.tracts.bed --output examples/sstar2.example.nean.match.rate.bed
sstar2 match --vcf examples/sstar2.example.biallelic.snps.vcf.gz --tgt examples/tgt.samples.list --src examples/den.samples.list --tract-file examples/sstar2.example.inferred.tracts.bed --output examples/sstar2.example.den.match.rate.bed
sstar2 assign --match-rate examples/sstar2.example.nean.match.rate.bed examples/sstar2.example.den.match.rate.bed --source-name nean den --output-prefix examples/sstar2.example
```
