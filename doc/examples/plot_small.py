import gaussdca
import pylab as plt

results = gaussdca.run('../../gaussdca/test/data/small.a3m')

plt.imshow(results['gdca_corr'])
plt.colorbar()
plt.show()

