from CADRE.CADRE_assembly import CADRE

from pylab import *
import cPickle

class generate_plots(object):
    
    def __init__(self, testAssembly, assembly):
        self.assembly = assembly
        self.npts = 6 # This value was specified in kwargs
        self.set_vals()
        self.assembly.run()
        #self.__initFigures()
        #self.__drawAxes()
        #self.plot('plot')
        #self.plot0('plot')
        #self.savePickled('data')

        self.myplot()
        self.plot_from_pkl()
        
    def set_vals(self):
        for comp in assembly.list_components():
            for input_name in assembly.get(comp).list_inputs():
                for var_name in testAssembly:
                    if var_name[0] == '5' or not var_name[0].isdigit():
                        if var_name[0].isdigit(): #do one at a time
                            var_name_short = var_name[2:]
                        if input_name == var_name_short:
                            try:
                                if isinstance(assembly.get(comp).get(input_name), float) and shape(testAssembly[var_name]) == (1,):
                                    assembly.get(comp).set(input_name, testAssembly[var_name][0])
                                elif isinstance(assembly.get(comp).get(input_name), ndarray) and shape(testAssembly[var_name]) == (1,):
                                    for i in range(len(assembly.get(comp).get(input_name))):
                                        assembly.get(comp).set(input_name, testAssembly[var_name]) #This may not work, but probably does
                                else:
                                    assembly.get(comp).set(input_name, testAssembly[var_name])
                                break
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
        time = linspace(0,43200.0,1500)
        xlim = [time.min()/3600.0, time.max()/3600.0]
        
        for curve in self.curves:
            f, i, j, r, c, scl, var, yLabel, ylim, lines = curve

            dx = 0 * (xlim[1] - xlim[0])
            dy = 0 * (ylim[1] - ylim[0])
            fig = figure(f)
            ax = fig.add_subplot(r, c, i*c + j + 1)
            if var[0].isdigit():
                var = var[2:]
            if var == 'gamma':
                var = 'Gamma'
            if len(assembly.get(assembly.varnames[var]).get(var).shape) is 1:
                curve[-1] = [plot(time/3600.0, time, 'k')[0]]
            else:
                n = assembly.get(assembly.varnames[var]).get(var).shape[0]
                curve[-1] = [plot(time/3600.0, time, 'k')[0] for k in range(n)]
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
        cPickle.dump(assembly.list_components(), f, 2)#cPickle.dump(self.main.val.arr, f, 2)
        f.close()
    
    def plot(self, filename):
        for curve in self.curves:
            f, i, j, r, c, scl, var, ylabel, ylim, lines = curve
            if var[0].isdigit():
                var_short = var[2:]
            if var_short == 'gamma':
                var_short = 'Gamma'
            fig = figure(f)
            ax = fig.add_subplot(r, c, i*c + j + 1)
            if len(assembly.get(assembly.varnames[var_short]).get(var_short).shape) is 1:
                lines[0].set_ydata(assembly.get(assembly.varnames[var_short]).get(var_short)*scl)
            else:
                n = assembly.get(assembly.varnames[var_short]).get(var_short).shape[0]
                for k in range(n):
                    lines[k].set_ydata(assembly.get(assembly.varnames[var_short]).get(var_short)[k,:]*scl)
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

            self.addSubplot(1,[180.0/pi*assembly.get(assembly.varnames['Gamma']).get('Gamma')], 'Roll angle [deg]', [0,90])
            self.addSubplot(4,[assembly.get(assembly.varnames['gain']).get('gain')], 'Gain', [0,2])
            self.addSubplot(7,[1e-3*assembly.get(assembly.varnames['Data']).get('Data')[0,:]], 'Total data [Gb]', [0,10])
            self.addSubplot(10,[assembly.get(assembly.varnames['CommLOS']).get('CommLOS')], 'Comm. line of sight', [0, 1])
            
            self.addSubplot(2,[assembly.get(assembly.varnames['Isetpt']).get('Isetpt')[k,:] for k in range(12)], 'Set pt. current [A]', [0,0.4])
            self.addSubplot(5,[assembly.get(assembly.varnames['temperature']).get('temperature')[k,:] for k in range(5)], 'Temp. [k]', [230,350])
            self.addSubplot(8,[assembly.get(assembly.varnames['V_sol']).get('V_sol')[k,:] for k in range(12)], 'Cell voltage [V]', [-2,6])
            self.addSubplot(11,[assembly.get(assembly.varnames['LOS']).get('LOS')], 'Solar line of sight', [0, 1])
            
            self.addSubplot(3,[assembly.get(assembly.varnames['P_comm']).get('P_comm')], 'Comm. power [W]', [-5,5])
            self.addSubplot(6,[assembly.get(assembly.varnames['P_sol']).get('P_sol')], 'Solar power [W]', [-8,8])
            self.addSubplot(9,[assembly.get(assembly.varnames['SOC']).get('SOC')[0,:]], 'State of charge', [0,1])
            self.addSubplot(12,[assembly.get(assembly.varnames['LOS']).get('LOS')], 'Solar line of sight', [0, 1])
            
            #pylab.savefig(filename + '-' + pre + '.png')
            savefig(filename + '-' + pre + '.pdf')

    def addSubplot(self, index, data, yLabel, ylim):
        time = linspace(0,43200.0,1500)

        xlim = [time.min()/3600.0, time.max()/3600.0]
        
        ax = subplot(self.r, self.c, index)
        for datum in data:
            plot(time/3600.0, datum, 'k')
        
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



    def myplot(self):
        time = linspace(0,43200.0,1500)
        xlim = [time.min()/3600.0, time.max()/3600.0]

        figure()
        
        subplot(4,3,1)
        xlabel = "Time"
        ylabel = "Roll Angle [deg]"
        ylim([0,90])
        plot(time, 180./pi*assembly.get(assembly.varnames['Gamma']).get('Gamma'))
        
        subplot(4,3,2)
        xlabel = "Time"
        ylabel = "Set pt. current [A]"
        ylim([0,0.4])
        for k in range (12):
            plot(time, assembly.get(assembly.varnames['Isetpt']).get('Isetpt')[k,:])
        
        subplot(4,3,3)
        xlabel = "Time"
        ylabel = "Comm. power [W]"
        ylim([-5,5])
        plot(time, assembly.get(assembly.varnames['P_comm']).get('P_comm'))

        subplot(4,3,4)
        xlabel = "Time"
        ylabel = "Gain"
        #ylim([0,2])
        plot(time, assembly.get(assembly.varnames['gain']).get('gain'))
        
        subplot(4,3,5)
        xlabel = "Time"
        ylabel = "Temp [K]"
        ylim([230,350])
        for k in range(5):
            plot(time, assembly.get(assembly.varnames['temperature']).get('temperature')[k])
        
        subplot(4,3,6)
        xlabel = "Time"
        ylabel = "Solar power [W]"
        ylim([-8,8])
        plot(time, assembly.get(assembly.varnames['P_sol']).get('P_sol'))
        
        subplot(4,3,7)
        xlabel = "Time"
        ylabel = "Total data [Gb]"
        #ylim([0,10])
        plot(time, 1e-3*assembly.get(assembly.varnames['Data']).get('Data')[0,:])
        
        subplot(4,3,8)
        xlabel = "Time"
        ylabel = "Cell voltage [V]"
        #ylim([-2,6])
        for k in range(12):
            plot(time, assembly.get(assembly.varnames['V_sol']).get('V_sol')[k,:])
        
        subplot(4,3,9)
        xlabel = "Time"
        ylabel = "State of charge"
        ylim([0,1])
        plot(time, assembly.get(assembly.varnames['SOC']).get('SOC')[0,:])
        
        subplot(4,3,10)
        xlabel = "Time"
        ylabel = "Comm. line of sight"
        #ylim([0,1])
        plot(time, assembly.get(assembly.varnames['CommLOS']).get('CommLOS'))
        
        subplot(4,3,11)
        xlabel = "Time"
        ylabel = "Solar line of sight"
        #ylim([0,1])
        plot(time, assembly.get(assembly.varnames['LOS']).get('LOS'))
        
        subplot(4,3,12)
        xlabel = "Time"
        ylabel = "Solar line of sight"
        #ylim([0,1])
        plot(time, assembly.get(assembly.varnames['LOS']).get('LOS'))
        
        savefig('myplot.pdf')


    def plot_from_pkl(self):
        time = linspace(0,43200.0,1500)
        xlim = [time.min()/3600.0, time.max()/3600.0]
        
        figure()
        
        subplot(4,3,1)
        xlabel = "Time"
        ylabel = "Roll Angle [deg]"
        ylim([0,90])
        plot(time, 180./pi*testAssembly['5:gamma'])
        
        subplot(4,3,2)
        xlabel = "Time"
        ylabel = "Set pt. current [A]"
        ylim([0,0.4])
        for k in range(12):
            plot(time, testAssembly['5:Isetpt'][k,:])
        
        subplot(4,3,3)
        xlabel = "Time"
        ylabel = "Comm. power [W]"
        ylim([-5,5])
        plot(time, testAssembly['5:P_comm'])
        
        subplot(4,3,4)
        xlabel = "Time"
        ylabel = "Gain"
        ylim([0,2])
        plot(time, testAssembly['5:gain'])
        
        subplot(4,3,5)
        xlabel = "Time"
        ylabel = "Temp [K]"
        ylim([230,350])
        for k in range(5):
            plot(time, testAssembly['5:temperature'][k])
        
        subplot(4,3,6)
        xlabel = "Time"
        ylabel = "Solar power [W]"
        ylim([-8,8])
        plot(time, testAssembly['5:P_sol'])
        
        subplot(4,3,7)
        xlabel = "Time"
        ylabel = "Total data [Gb]"
        ylim([0,10])
        plot(time, 1e-3*testAssembly['5:Data'][0,:])
        
        subplot(4,3,8)
        xlabel = "Time"
        ylabel = "Cell voltage [V]"
        ylim([-2,6])
        for k in range(12):
            plot(time, testAssembly['5:V_sol'][k,:])
        
        subplot(4,3,9)
        xlabel = "Time"
        ylabel = "State of charge"
        ylim([0,1])
        plot(time, testAssembly['5:SOC'][0,:])
        
        subplot(4,3,10)
        xlabel = "Time"
        ylabel = "Comm. line of sight"
        ylim([0,1])
        plot(time, testAssembly['5:CommLOS'])
        
        subplot(4,3,11)
        xlabel = "Time"
        ylabel = "Solar line of sight"
        ylim([0,1])
        plot(time, testAssembly['5:LOS'])
        
        subplot(4,3,12)
        xlabel = "Time"
        ylabel = "Solar line of sight"
        ylim([0,1])
        plot(time, testAssembly['5:LOS'])
        
        savefig('plot_from_pkl.pdf')







if __name__ == "__main__":
    testAssembly = cPickle.load(open('data1346.pkl', 'r'))
    assembly = CADRE()
    generate_plots(testAssembly, assembly)