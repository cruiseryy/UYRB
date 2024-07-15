# adopted from Julliane Quinn's code 
# https://github.com/julianneq/RedRiver_RivalFramings/blob/master/PaperFigures/makeFigure4.py
import numpy as np
import matplotlib
from matplotlib import pyplot as plt
from matplotlib import patheffects as pe

class parallel_coords:
    def __init__(self,
                 mins,
                 maxs,
                 xlabels,
                 precision):
        self.mins = mins
        self.maxs = maxs
        self.xlabels = xlabels
        self.precision = precision
        return
    
    def plot_(self, ax, data, cmap, shading):
        toplabels = []
        botlabels = []
        for i in range(len(self.xlabels)):
            if self.precision[i] != 0:
                toplabels.append(str(np.round(self.maxs[i], self.precision[i])))
                botlabels.append(str(np.round(self.mins[i], self.precision[i])))
            else:
                toplabels.append(str(int(self.maxs[i])))
                botlabels.append(str(int(self.mins[i])))
            if self.maxs[i] < 0:
                toplabels[i] = toplabels[i][1:]
            if self.mins[i] < 0:
                botlabels[i] = botlabels[i][1:]
       
        cmap = matplotlib.cm.get_cmap(cmap)
        scaled = np.zeros(data.shape)
        for j in range(data.shape[1]):
            scaled[:, j] = (data[:, j] - self.mins[j]) / (self.maxs[j] - self.mins[j])
        
        # his_obj = np.array([-324.4417, -5.8861, 0.0788, -0.5763])
        # his_obj_scaled = (his_obj - mins)/(maxs - mins)
        # dev = np.zeros([data.shape[0], ])
        # for i in range(data.shape[0]):
        #     dev[i] = np.sum((scaled[i, :] - his_obj_scaled)**2)
        # dev = 1 - (dev - np.min(dev)) / (np.max(dev) - np.min(dev))

        xs = np.arange(data.shape[1])
        for i in range(data.shape[0]):
            ys = scaled[i, :]
            ax.plot(xs, ys, c=cmap(0.2+0.7*shading[i]), linewidth=2, alpha = 0.25)

        # min_dev_idx = np.argmax(dev)
        # ax.plot(xs, scaled[min_dev_idx, :], c=cmap(0.9), linewidth=4)
        # ax.plot(xs, his_obj_scaled, c='black', linewidth=4)
        for i in range(data.shape[1]):
            ax.axvline(i, linewidth=1.5, color='k', zorder=10)

        ax.set_xticks(xs)
        ax.set_xlim([xs[0]-0.15, xs[-1]+0.15])
        ax.set_ylim([0, 1])
        
        ax.set_yticks([])
        ax.set_xticks([])
        ax.set_xticklabels([])

        ax.spines["top"].set_visible(False)
        ax.spines["bottom"].set_visible(False)
        ax.spines["left"].set_visible(False)
        ax.spines["right"].set_visible(False)

        for i, x in enumerate(xs):
            ax.text(x, 1.02, toplabels[i], ha='center', va='bottom', fontsize=16)
            ax.text(x, -0.02, botlabels[i], ha='center', va='top', fontsize=16)


        return 
