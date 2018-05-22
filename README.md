# pyGaussDCA

Python re-implementation of GaussDCA by ["Fast and accurate multivariate
Gaussian modeling of protein families: Predicting residue contacts and
protein-interaction partners"][paper]
by Carlo Baldassi, Marco Zamparo, Christoph Feinauer, Andrea Procaccini,
Riccardo Zecchina, Martin Weigt and Andrea Pagnani, (2014)
PLoS ONE 9(3): e92721. doi:10.1371/journal.pone.0092721.
It is based on the [original code][original] in Julia


[paper]: http://www.plosone.org/article/info%3Adoi%2F10.1371%2Fjournal.pone.0092721
[original]: https://github.com/carlobaldassi/GaussDCA.jl

## What is different?
It is written as a Python library to be more easily integrated with our software. 

## How do I install it?
First install [Pythran][pythran] 0.8.5 or higher. Then, just use pip to install pyGaussDCA.

[pythran]:https://github.com/serge-sans-paille/pythran

## How do I use it?

Pass the path to an alignment file to the function. The supported formats are A3M, FASTA, and ALN, without line wrap: every sequence spans one single line.

```
import gaussdca

results = gaussdca.run('/path/to/a3m')
```

## TODO:
  * Example and doc
  * Tests
  * CLI util
  * Pypi
  
