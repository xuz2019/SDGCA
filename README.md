Code of "Similarity and Dissimilarity Guided Co-association Matrix Construction for Ensemble Clustering".

The hyperparameters $\lambda$, $\eta$, $\theta$ can be tuned for better performance. In our experiment, the hyperparameters for different datasets are set to

```
lambda = ('Ecoli', 0.09, 'GLIOMA', 0.02, 'Aggregation', 0.08, 'MF', 0.05, 'IS', 0.03, 'MNIST', 0.07, 'Texture',0.04, 'SPF',0.06, 'ODR', 0.06, 'LS', 0.18, 'ISOLET', 0.06, 'USPS', 0.06, 'orlraws10P', 0.95, 'BBC', 0.06, 'lung', 0.75);
```
```
eta =    ('Ecoli', 0.65, 'GLIOMA', 0.8,  'Aggregation', 0.65, 'MF', 0.75, 'IS', 0.9,  'MNIST', 0.95, 'Texture', 1,   'SPF', 0.7, 'ODR', 0.95, 'LS', 0.7,  'ISOLET', 0.95, 'USPS', 0.95, 'orlraws10P', 1.01, 'BBC', 1.01, 'lung', 0.75);
```
```
theta  = ('Ecoli', 0.75, 'GLIOMA', 0.6,  'Aggregation', 0.7,  'MF', 0.95, 'IS', 0.95, 'MNIST', 0.95, 'Texture', 1,   'SPF', 0.9, 'ODR', 0.95, 'LS', 0.6,  'ISOLET', 0.95, 'USPS', 0.95, 'orlraws10P', 1.01, 'BBC', 1.01, 'lung', 0.6);
```
