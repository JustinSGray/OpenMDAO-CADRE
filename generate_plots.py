from CADRE.CADRE_assembly import CADRE

from pylab import *
import cPickle

class generate_plots(object):
    
    def __init__(self, testAssembly, assembly):
        self.assembly = assembly
        self.npts = 6 # This value was specified in kwargs
        self.set_vals()
        self.__initFigures()
        self.__drawAxes()
        self.plot('plot')
        self.plot0('plot')
        self.savePickled('data')
    
    def set_vals(self):
        for comp in assembly.list_components():
            for input_name in assembly.get(comp).list_inputs():
                for var_name in testAssembly:
                    if var_name[0].isdigit():
                        var_name_short = var_name[2:]
                    if input_name == var_name_short:
                        try:
                            if shape(testAssembly[var_name]) == (1,):
                                assembly.get(comp).set(input_name, testAssembly[var_name][0])
                            else:
                                assembly.get(comp).set(input_name, testAssembly[var_name])
                        except RuntimeError:
                            pass
               
    def __initFigures(self):
        fig_width = 20  # width in inches
        fig_height = 12  # height in inches
        fig_size =  [fig_width,fig_height]
        params = {#'backend': 'Agg',
            'axes.labelsize': 16,
            'axes.titlesize': 16,
            'font.size': 16,
            'xtick.labelsize': 16,
            'ytick.labelsize': 16,
            'figure.figsize': fig_size,
            'savefig.dpi' : 600,
            'font.family': 'sans-serif',
            'axes.linewidth' : 0.5,
            'xtick.major.size' : 4,
            'ytick.major.size' : 4,
            'font.size' : 16
        }
        rcParams.update(params)
        
        for k in range(self.npts+1):
            fig = figure(k)
            fig.subplots_adjust(wspace=0.6, hspace=0.6)

        curves = []
        for f in [0]:
            r,c = 3,6
            for pt in range(self.npts):
                pre = str(pt)
                curves.append([f, 0, pt, r, c, 1e-3, pre+':Data', 'Total data [Gb]', [0,10],None])
                curves.append([f, 1, pt, r, c, 1.0, pre+':V_sol', 'Cell voltage [V]', [-2,8],None])
                curves.append([f, 2, pt, r, c, 1.0, pre+':SOC', 'State of charge', [0,1],None])
        for f in range(1,self.npts+1):
            r,c = 4,3
            pre = str(f-1)
            curves.append([f, 0, 0, r, c, 180.0/pi, pre+':gamma', 'Roll angle [deg]', [0,90],None])
            curves.append([f, 0, 1, r, c, 1.0, pre+':Isetpt', 'Set pt. current [A]', [0,0.4],None])
            curves.append([f, 0, 2, r, c, 1.0, pre+':P_comm', 'Comm. power [W]', [-5,5],None])
            curves.append([f, 1, 0, r, c, 1.0, pre+':gain', 'Gain', [0,2],None])
            curves.append([f, 1, 1, r, c, 1.0, pre+':temperature', 'Temp. [k]', [230,350],None])
            curves.append([f, 1, 2, r, c, 1.0, pre+':P_sol', 'Solar power [W]', [-8,8],None])
            curves.append([f, 2, 0, r, c, 1e-3, pre+':Data', 'Total data [Gb]', [0,10],None])
            curves.append([f, 2, 1, r, c, 1.0, pre+':V_sol', 'Cell voltage [V]', [-2,8],None])
            curves.append([f, 2, 2, r, c, 1.0, pre+':SOC', 'State of charge', [0,1],None])
            curves.append([f, 3, 0, r, c, 1.0, pre+':CommLOS', 'Comm. line of sight', [0, 1],None])
            curves.append([f, 3, 1, r, c, 1.0, pre+':LOS', 'Solar line of sight', [0, 1],None])
            curves.append([f, 3, 2, r, c, 1.0, pre+':LOS', 'Solar line of sight', [0, 1],None])
                
        self.curves = curves
    
    def __drawAxes(self):
        t = linspace(0,43200.,1500) # THIS IS JUST A GUESS BASED ON OUTPUT.TXT
        xlim = [t.min()/3600.0, t.max()/3600.0]
        
        for curve in self.curves:
            f, i, j, r, c, scl, var, yLabel, ylim, lines = curve

            dx = 0 * (xlim[1] - xlim[0])
            dy = 0 * (ylim[1] - ylim[0])
            fig = figure(f)
            ax = fig.add_subplot(r, c, i*c + j + 1)
            if len(var) is 1: # used to be: if len(self.main.val(var).shape) is 1:
                curve[-1] = [plot(t/3600.0, t, 'k')[0]]
            else:
                n = len(var[0])# used to be: n = self.main.val(var).shape[0]
                curve[-1] = [plot(t/3600.0, t, 'k')[0] for k in range(n)]
            ax.set_xlim((xlim[0]-dx,xlim[1]+dx))
            ax.set_ylim((ylim[0]-dy,ylim[1]+dy))
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
            ax.spines['bottom'].set_position(('outward', 15))
            ax.spines['left'].set_position(('outward', 18))
            ax.yaxis.set_ticks_position('left')
            ax.xaxis.set_ticks_position('bottom')

            xlabel('time (hr)')
            ylabel(yLabel)
            
            xticks(linspace(xlim[0],xlim[1],5), size=12)
            yticks(linspace(ylim[0],ylim[1],5), size=12)
    
    def savePickled(self, filename):
        f = open(filename+'.pkl','wb')
        cPickle.dump(assembly, f, 2)#cPickle.dump(self.main.val.arr, f, 2)
        f.close()
    
    def plot(self, filename):
        for curve in self.curves:
            f, i, j, r, c, scl, var, ylabel, ylim, lines = curve
            if var[0].isdigit(): #THIS MAY NOT WORK............................................................
                var = var[2:]
            fig = figure(f)
            ax = fig.add_subplot(r, c, i*c + j + 1)
            if len(assembly.get(var).shape) is 1:#if len(self.main.val(var).shape) is 1:
                lines[0].set_ydata(assembly.get(var)*scl)#lines[0].set_ydata(self.main.val(var)*scl)
            else:
                n = assembly.get(var).shape[0]#n = self.main.val(var).shape[0]
                for k in range(n):
                    #lines = assembly.get(var)*scl#lines[k].set_ydata(self.main.val(var)[k,:]*scl)
                    lines[k].set_ydata(assembly.get(var)[k,:])
        for f in [0]:
            fig = figure(f)
            fig.canvas.draw()
            #savefig(filename + '.png')
            savefig(filename + '.pdf')
        for f in range(1,self.npts+1):
            pre = str(f-1)
            fig = figure(f)
            fig.canvas.draw()
            #savefig(filename + '-' + pre + '.png')
            savefig(filename + '-' + pre + '.pdf')

    def plot0(self, filename):        

        fig = figure(0)
        
        curves = []
        
        
        
        for pt in range(self.npts):
            figure(pt+1)
            self.r, self.c = 4, 3
            pre = str(pt)
            
            self.addSubplot(1,[180.0/pi*assembly.get(pre+':gamma')], 'Roll angle [deg]', [0,90])
            self.addSubplot(4,[assembly.get(pre+':gain')], 'Gain', [0,2])
            self.addSubplot(7,[1e-3*assembly.get(pre+':Data')[0,:]], 'Total data [Gb]', [0,10])
            self.addSubplot(10,[assembly.get(pre+':CommLOS')], 'Comm. line of sight', [0, 1])
            
            self.addSubplot(2,[assembly.get(pre+':Isetpt')[k,:] for k in range(12)], 'Set pt. current [A]', [0,0.4])
            self.addSubplot(5,[assembly.get(pre+':temperature')[k,:] for k in range(5)], 'Temp. [k]', [230,350])
            self.addSubplot(8,[assembly.get(pre+':V_sol')[k,:] for k in range(12)], 'Cell voltage [V]', [-2,6])
            self.addSubplot(11,[assembly.get(pre+':LOS')], 'Solar line of sight', [0, 1])
            
            self.addSubplot(3,[assembly.get(pre+':P_comm')], 'Comm. power [W]', [-5,5])
            self.addSubplot(6,[assembly.get(pre+':P_sol')], 'Solar power [W]', [-8,8])
            self.addSubplot(9,[assembly.get(pre+':SOC')[0,:]], 'State of charge', [0,1])
            self.addSubplot(12,[assembly.get(pre+':LOS')], 'Solar line of sight', [0, 1])
            
            #pylab.savefig(filename + '-' + pre + '.png')
            savefig(filename + '-' + pre + '.pdf')

    def addSubplot(self, index, data, yLabel, ylim):
        t = linspace(0,43200.,1500) # THIS IS JUST A GUESS BASED ON OUTPUT.TXT

        xlim = [t.min()/3600.0, t.max()/3600.0]
        
        ax = subplot(self.r, self.c, index)
        for datum in data:
            plot(t/3600.0, datum, 'k')
        
        ax.set_xlim((xlim[0],xlim[1]))
        ax.set_ylim((ylim[0],ylim[1]))
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['bottom'].set_position(('outward', 15))
        ax.spines['left'].set_position(('outward', 18))
        ax.yaxis.set_ticks_position('left')
        ax.xaxis.set_ticks_position('bottom')
        
        xlabel('time (hr)')
        ylabel(yLabel)
        
        xticks(linspace(xlim[0],xlim[1],5), size=12)
        yticks(linspace(ylim[0],ylim[1],5), size=12)

if __name__ == "__main__":
    testAssembly = cPickle.load(open('data1346.pkl', 'r'))
    assembly = CADRE()
    generate_plots(testAssembly, assembly)