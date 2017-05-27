import sys
import pandas as pd
from matplotlib import pyplot as plt

fil = sys.argv[1]

data = pd.read_csv(fil, sep = '\t', header= None, dtype = {0: str})

pivoted = data.pivot_table(index = [0,1], columns = 2, values = 3)

pivoted.columns.name = ''
pivoted.index.name = 'Genome numbers'

pivoted.plot(y = pivoted.columns, kind = 'bar', use_index=True, colormap = 'Pastel1', figsize = (20,22))

plt.savefig(fil[:-4]+'.png')
