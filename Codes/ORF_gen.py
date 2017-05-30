import sys
import pandas as pd
from matplotlib import pyplot as plt

fil = sys.argv[1]

spec_dic = {'05': 'Chlamydia trachomatis', '08': 'Dictyoglomus turgidum', '14': 'Mycobacterium tuberculosis', '16': 'Rhodopirellula baltica', '25': 'Saccharomyces cerevisiae'}

data = pd.read_csv(fil, sep = '\t', header= None, dtype = {0: str})

data.replace({0: spec_dic}, inplace = True)

pivoted = data.pivot_table(index = [0,1], columns = 2, values = 3)

pivoted.columns.name = ''
pivoted.index.name = 'Genome numbers'

pivoted.plot(y = pivoted.columns, kind = 'bar', use_index=True, colormap = 'Pastel1', figsize = (10,10), legend = False)

plt.tight_layout()
#plt.legend(loc='lower left', bbox_to_anchor=(-0.05, -0.3))
plt.xlabel('Species')
plt.savefig(fil[:-4]+'.png')
