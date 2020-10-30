#to Do's:   - get rid of functions_single.py
#           - make a simulation of a larger protein and check if resid change works

#Imports

import MDAnalysis as mda
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from .functions_single import prepare_dihedrals, angle_to_list, angle_trans, angle_diff
from .functions_capped import prepare_capped_dihedrals
from .functions_just_top import prepare_top_dihedrals
from .functions_multi import prepare_all_dihedrals
import click
import seaborn as sns
import statistics
import numpy as np
import scipy.stats as st
import matplotlib.animation as anim

#Make it a little more artsy

#Set up the command line options
#Only topology and trajectories are needed

@click.group()
@click.option('--top', help='path to topology.gro or topology.pdb')
@click.option('--trj', help='path to trajectories.xtc or trajectories.dcd', default='')
@click.option('--res', help='residue number (default = 2)', default=2)
@click.option('--id', help='index of residue number in .pdb file (default = 4)', default=4)
@click.option('--ma', help='animates all residues in a protein', is_flag=True)
@click.option('--start', help='Only works when --ma flag is set, gives the number of the first residues that should be animated (default = 2)', type=int, required=False, default=2)
@click.option('--end', help='Only works when --ma flag is set, gives the number of the last residues that should be animated (default = 2)', type=int, required=False, default=2)
@click.pass_context

#Calculates the angles from the input topology and trajectories and store it in ctx.obj
#Find more about angle calculations in functions.py





def main(ctx, top: str, trj: str, res: str, id:str, ma:bool, start:int, end:int):
    if trj == '':
        ctx.obj['trj'] = 'empty'
        phi_gen, psi_gen, phi_pro, psi_pro, phi_gly, psi_gly, phi_pre, psi_pre = prepare_top_dihedrals(top)
        ctx.obj['psi_gen'] = psi_gen
        ctx.obj['phi_gen'] = phi_gen
        ctx.obj['psi_pro'] = psi_pro
        ctx.obj['phi_pro'] = phi_pro
        ctx.obj['psi_gly'] = psi_gly
        ctx.obj['phi_gly'] = phi_gly
        ctx.obj['psi_pre'] = psi_pre
        ctx.obj['phi_pre'] = phi_pre
    else:
        ctx.obj['trj'] = 'filled'
        u = mda.Universe(top, trj)
        if ma == True:
            if end < start:
                print('Wait a minute, that is not possible!')
            else:
                phi, psi, res_name_list = prepare_all_dihedrals(top, u, res, id, start, end+1)
                ctx.obj['psi'] = psi
                ctx.obj['phi'] = phi
                ctx.obj['res'] = res_name_list
        elif 'ACE' in str(u.residues[0]):
            phi, psi = prepare_capped_dihedrals(top, u, res, id)
            ctx.obj['psi'] = psi
            ctx.obj['phi'] = phi
        #This is only for uncapped, single amino acids
        #Only .gro/.xtc files are working here
        else:
            phi, psi = prepare_dihedrals(top, u)
            ctx.obj['psi'] = psi
            ctx.obj['phi'] = phi



#Set up the commands
#Plot the data and give some informative visualizations

@main.command()
@click.pass_context
def plot(ctx):
    # Sets up the figure and subplots
    if ctx.obj['trj'] == 'empty':

        im_gen = plt.imread('/home/josh/PycharmProjects/saarama_project/saarama/bg/rama_bg.png')
        im_gly = plt.imread('/home/josh/PycharmProjects/saarama_project/saarama/bg/gly_bg.png')
        im_pro = plt.imread('/home/josh/PycharmProjects/saarama_project/saarama/bg/pro_bg.png')
        im_pre = plt.imread('/home/josh/PycharmProjects/saarama_project/saarama/bg/pre_bg.png')

        fig = plt.figure(constrained_layout=True)
        gs = gridspec.GridSpec(2, 2, figure=fig)

        ax1 = fig.add_subplot(gs[0, 0])
        ax1.imshow(im_gen, extent=(-180,180,-180,180))
        ax1.scatter(ctx.obj['phi_gen'], ctx.obj['psi_gen'], s=15, color='dimgray')
        ax1.plot([0, 0], [-180, 180], c='k', alpha=0.3)
        ax1.plot([-180, 180], [0, 0], c='k', alpha=0.3)
        ax1.set_xlim(-180, 180)
        ax1.set_ylim(-180, 180)
        ax1.set_xlabel('φ')
        ax1.set_ylabel('ψ')
        ax1.set_title('Ramachandran plot')


        ax2 = fig.add_subplot(gs[0, 1])
        ax2.scatter(ctx.obj['phi_gly'], ctx.obj['psi_gly'], s=15, color='dimgray')
        ax2.imshow(im_gly, extent=(-180,180,-180,180))
        ax2.plot([0, 0], [-180, 180], c='k', alpha=0.3)
        ax2.plot([-180, 180], [0, 0], c='k', alpha=0.3)
        ax2.set_xlim(-180, 180)
        ax2.set_ylim(-180, 180)
        ax2.set_xlabel('φ')
        ax2.set_ylabel('ψ')
        ax2.set_title('Ramachandran plot of Glycine')

        ax3 = fig.add_subplot(gs[1, 0])
        ax3.scatter(ctx.obj['phi_pro'], ctx.obj['psi_pro'], s=15, color='dimgray')
        ax3.imshow(im_pro, extent=(-180,180,-180,180))
        ax3.plot([0, 0], [-180, 180], c='k', alpha=0.3)
        ax3.plot([-180, 180], [0, 0], c='k', alpha=0.3)
        ax3.set_xlim(-180, 180)
        ax3.set_ylim(-180, 180)
        ax3.set_xlabel('φ')
        ax3.set_ylabel('ψ')
        ax3.set_title('Ramachandran plot of Proline')

        ax4 = fig.add_subplot(gs[1, 1])
        ax4.scatter(ctx.obj['phi_pre'], ctx.obj['psi_pre'], s=15, color='dimgray')
        ax4.imshow(im_pre, extent=(-180,180,-180,180))
        ax4.plot([0, 0], [-180, 180], c='k', alpha=0.3)
        ax4.plot([-180, 180], [0, 0], c='k', alpha=0.3)
        ax4.set_xlim(-180, 180)
        ax4.set_ylim(-180, 180)
        ax4.set_xlabel('φ')
        ax4.set_ylabel('ψ')
        ax4.set_title('Ramachandran plot of Pre-proline')

        plt.show()

    elif ctx.obj['trj'] == 'filled':

        plt.style.use('seaborn-darkgrid')

        fig = plt.figure(constrained_layout=True)
        gs = gridspec.GridSpec(3, 2, figure=fig)
        fig.suptitle('Number of angles: ' + str(len(ctx.obj['psi'])), fontsize=12)

        # Scatter plot

        ax1 = fig.add_subplot(gs[0, 0])
        ax1.scatter(ctx.obj['phi'], ctx.obj['psi'], s=15, color='dimgray')
        ax1.plot([0, 0], [-180, 180], c='k', alpha=0.3)
        ax1.plot([-180, 180], [0, 0], c='k', alpha=0.3)
        ax1.set_xlim(-180, 180)
        ax1.set_ylim(-180, 180)
        ax1.set_xlabel('φ')
        ax1.set_ylabel('ψ')
        ax1.set_title('Ramachandran plot of a single amino acid')

        # Contour plot

        ax2 = fig.add_subplot(gs[0, 1])
        sns.kdeplot(ctx.obj['phi'], ctx.obj['psi'], ax=ax2, cmap='Reds', shade=True, shade_lowest=False)
        sns.rugplot(ctx.obj['phi'], color='k', ax=ax2)
        sns.rugplot(ctx.obj['psi'], color='k', vertical=True, ax=ax2)
        ax2.plot([0, 0], [-180, 180], c='k', alpha=0.3)
        ax2.plot([-180, 180], [0, 0], c='k', alpha=0.3)
        ax2.set_xlabel('φ')
        ax2.set_ylabel('ψ')
        ax2.set_title('Contour plot of a single amino acid')

        '''
        N_bins = 120

        counts, xedges, yedges, im = ax2.hist2d(ctx.obj['phi'], ctx.obj['psi'], bins=N_bins, density=True, cmap='plasma')
        fig.colorbar(im, ax=ax2)
        '''
        # Calculates angle differences

        psi_trans = angle_trans(ctx.obj['psi'])
        phi_trans = angle_trans(ctx.obj['phi'])
        psi_diff = angle_diff(psi_trans)
        phi_diff = angle_diff(phi_trans)

        combined_list = psi_diff + phi_diff

        # Angle difference over time for Phi

        ax3 = fig.add_subplot(gs[1, 0])
        # average_phi = sum(ctx.obj['phi'])/len(ctx.obj['phi'])
        # ax3.plot([0, len(ctx.obj['phi'])], [average_phi, average_phi], color='k', alpha=0.75, label='φ-avg: '+str(round(average_phi, 2)))
        # ax3.plot(range(0, len(ctx.obj['phi'])), ctx.obj['phi'], color='indianred')
        average_phi = sum(phi_diff) / len(phi_diff)
        stdev_phi = statistics.stdev(phi_diff)
        ax3.plot([0, len(phi_diff)], [average_phi, average_phi], color='k', alpha=0.75,
                 label='φ-avg: ' + str(round(average_phi, 2)) + ', stdev: ' + str(round(stdev_phi, 2)))
        ax3.plot(phi_diff, color='indianred', alpha=0.5)
        ax3.scatter(range(len(phi_diff)), phi_diff, marker='x', s=5, color='k')
        ax3.set_ylabel('Angel difference')
        ax3.set_xlabel('Timeframe')
        ax3.set_title('φ-difference over time')
        ax3.set_ylim(min(combined_list) - 5, max(combined_list) + 5)
        ax3.set_xlim(0, len(ctx.obj['phi']))
        plt.legend(loc='upper right', borderaxespad=0.)

        # Angle difference over time for Psi

        ax4 = fig.add_subplot(gs[1, 1])
        # average_psi = sum(ctx.obj['psi'])/len(ctx.obj['psi'])
        # ax4.plot([0, len(ctx.obj['psi'])], [average_psi, average_psi], color='k', alpha=0.75, label='ψ-avg: '+str(round(average_psi, 2)))
        # ax4.plot(range(0, len(ctx.obj['psi'])), ctx.obj['psi'], color='darkkhaki')
        average_psi = sum(psi_diff) / len(psi_diff)
        stdev_psi = statistics.stdev(psi_diff)
        ax4.plot([0, len(psi_diff)], [average_psi, average_psi], color='k', alpha=0.75,
                 label='ψ-avg: ' + str(round(average_psi, 2)) + ', stdev: ' + str(round(stdev_psi, 2)))
        ax4.plot(psi_diff, color='darkkhaki', alpha=0.5)
        ax4.scatter(range(len(psi_diff)), psi_diff, marker='x', s=5, color='k')
        ax4.set_ylabel('Angel difference')
        ax4.set_xlabel('Timeframe')
        ax4.set_title('ψ-difference over time')
        ax4.set_ylim(min(combined_list) - 5, max(combined_list) + 5)
        ax4.set_xlim(0, len(ctx.obj['phi']))
        plt.legend(loc='upper right', borderaxespad=0.)

        # 3D Density plot

        # bins_list_phi = list(range(int(min(psi_trans)), int(max(psi_trans)), 1))
        # bins_list_psi = list(range(int(min(phi_trans)), int(max(phi_trans)), 1))
        bins_list_phi = list(range(int(min(ctx.obj['phi'])), int(max(ctx.obj['phi'])), 1))
        bins_list_psi = list(range(int(min(ctx.obj['psi'])), int(max(ctx.obj['psi'])), 1))

        ax5 = fig.add_subplot(gs[2, 0], projection='3d')

        x = np.asarray(ctx.obj['phi'])
        y = np.asarray(ctx.obj['psi'])

        deltaX = (max(x) - min(x)) / 10
        deltaY = (max(y) - min(y)) / 10

        xmin = min(x) - deltaX
        xmax = max(x) + deltaX
        ymin = min(y) - deltaY
        ymax = max(y) + deltaY

        xx, yy = np.mgrid[xmin:xmax:100j, ymin:ymax:100j]

        positions = np.vstack([xx.ravel(), yy.ravel()])
        values = np.vstack([x, y])
        kernel = st.gaussian_kde(values)
        f = np.reshape(kernel(positions).T, xx.shape)

        # ax5.plot_wireframe(xx, yy, f, alpha=0.8)
        ax5.plot_surface(xx, yy, f, rstride=1, cstride=1, edgecolor='none', cmap='plasma')
        ax5.set_xlim(-180, 180)
        ax5.set_ylim(-180, 180)
        ax5.set_xlabel('φ')
        ax5.set_ylabel('ψ')
        ax5.set_zlabel('Density')
        ax5.set_title('Surface plot of angle distributions')
        ax5.view_init(20, 280)

        # Histogram/Density plot for angle distribution

        ax6 = fig.add_subplot(gs[2, 1])
        # sns.distplot(phi_trans, ax=ax5, bins=bins_list_phi, color = 'indianred', label='φ')
        # sns.distplot(psi_trans, ax=ax5, bins=bins_list_psi, color = 'darkkhaki', label='ψ')
        # ax5.set_xlim(0,360)
        # ax5.hist(ctx.obj['phi'], bins=bins_list_phi, alpha=0.75, color = 'indianred', label='φ')
        # ax5.hist(ctx.obj['psi'], bins=bins_list_psi, alpha=0.75, color = 'darkkhaki', label='ψ')
        sns.distplot(ctx.obj['phi'], ax=ax6, bins=bins_list_phi, color='indianred', label='φ')
        sns.distplot(ctx.obj['psi'], ax=ax6, bins=bins_list_psi, color='darkkhaki', label='ψ')
        ax6.set_xlim(-180, 180)
        ax6.set_title('Histogram/Density plot of angle distribution')
        plt.legend()

        '''
        #Polar plots that are not included yet

        ax6 = fig.add_subplot(gs[3, 0], projection='polar')
        bin_size = 20
        a, b = np.histogram(phi_trans, bins=np.arange(0, 360 + bin_size, bin_size))
        centers = np.deg2rad(np.ediff1d(b) // 2 + b[:-1])
        ax6.bar(centers, a, width=np.deg2rad(bin_size), bottom=0.0, color='.8', edgecolor='k')
        ax6.set_theta_zero_location("N")
        ax6.set_theta_direction(-1)

        ax7 = fig.add_subplot(gs[3, 1], projection='polar')
        bin_size = 20
        a, b = np.histogram(psi_trans, bins=np.arange(0, 360 + bin_size, bin_size))
        centers = np.deg2rad(np.ediff1d(b) // 2 + b[:-1])
        ax7.bar(centers, a, width=np.deg2rad(bin_size), bottom=0.0, color='.8', edgecolor='k')
        ax7.set_theta_zero_location("N")
        ax7.set_theta_direction(-1)
        '''

        plt.show()






#Print the data to stdout

@main.command()
@click.pass_context
def terminal(ctx):

    if ctx.obj['trj'] == 'empty':
        print('General torsion angles:')
        print('Phi: ', ctx.obj['phi_gen'])
        print('Psi: ', ctx.obj['psi_gen'])
        print('Glycine:')
        print('Phi: ', ctx.obj['phi_gly'])
        print('Psi: ', ctx.obj['psi_gly'])
        print('Proline:')
        print('Phi: ', ctx.obj['phi_pro'])
        print('Psi: ', ctx.obj['psi_pro'])
        print('Pre-proline:')
        print('Phi: ', ctx.obj['phi_pre'])
        print('Psi: ', ctx.obj['psi_pre'])
    elif ctx.obj['trj'] == 'filled':
        print('Torsion angles:')
        print('Phi: ', ctx.obj['phi'])
        print('Psi: ', ctx.obj['psi'])



@main.command()
@click.pass_context
def animate(ctx):

    phi = []
    psi = []

    if ctx.obj['trj'] == 'empty':
        phi = list(ctx.obj['phi_gen'])
        psi = list(ctx.obj['psi_gen'])
    elif ctx.obj['trj'] == 'filled':
        phi = ctx.obj['phi']
        psi = ctx.obj['psi']


    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    ax.set_xlim(-180, 180)
    ax.set_ylim(-180, 180)

    frames = len(phi)

    def animate(i):
        step = int(frames*0.01)
        i = int(i*step)
        j = int(i+step)
        x = phi[i:j]
        y = psi[i:j]
        ax.scatter(x, y, color='k', s=15)

    a = anim.FuncAnimation(fig, animate, frames=frames, interval=1, repeat=True, blit=False)
    ax.plot([0, 0], [-180, 180], c='k', alpha=0.2)
    ax.plot([-180, 180], [0, 0], c='k', alpha=0.2)
    plt.xlim(-180, 180)
    plt.ylim(-180, 180)
    plt.ylabel('ψ')
    plt.xlabel('φ')
    plt.show()

@main.command()
@click.pass_context
def multi_animate(ctx):

    phi = []
    psi = []

    phi = ctx.obj['phi']
    psi = ctx.obj['psi']
    res_name_list = ctx.obj['res']

    fig = plt.figure()

    columns = 0
    rows = 0

    if len(phi) < 4:
        rows = len(phi)
    else:
        rows = 4


    if len(phi) <= 4:
        columns = 1
    elif (len(phi) / 4) > 1:
        columns = int((len(phi) / 4))+1

    gs = gridspec.GridSpec(columns, rows)
    plot_names = ["ax" + str(i) for i in range(1, len(phi)+1)]
    axs = []

    i = 0

    for c, num in zip(plot_names, range(1, len(phi)+1)):
        axs.append(fig.add_subplot(gs[num - 1]))
        axs[-1].set_xlim(-180, 180)
        axs[-1].set_ylim(-180, 180)
        axs[-1].plot([0, 0], [-180, 180], c='k', alpha=0.2)
        axs[-1].plot([-180, 180], [0, 0], c='k', alpha=0.2)
        axs[-1].set_ylabel('ψ')
        axs[-1].set_xlabel('φ')
        axs[-1].set_title(res_name_list[i])
        i += 1

    def animate(i):
        step = int(frames*0.01)
        i = int(i*step)
        j = int(i+step)
        for k in range(0, len(phi)):
            x = phi[k][i:j]
            y = psi[k][i:j]
            axs[k].scatter(x,y,color='k', s=15)

    frames = len(phi[0])
    a = anim.FuncAnimation(fig, animate, frames=frames, interval=1, repeat=True, blit=False)

    #plt.tight_layout()
    plt.show()

def start():
    main(obj={})


if __name__ == '__main__':
    start()
