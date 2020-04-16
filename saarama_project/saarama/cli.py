#Imports

import MDAnalysis as mda
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from .functions import prepare_dihedrals, angle_to_list, angle_trans, angle_diff
import click
import seaborn as sns
import statistics

#Make it a little more artsy

plt.style.use('seaborn-darkgrid')

#Set up the command line options
#Only topology and trajectories are needed

@click.group()
@click.option('--top', help='path to topology.gro or topology.pdb')
@click.option('--trj', help='path to trajectories.xtc or trajectories.dcd')
@click.pass_context

#Calculates the angles from the input topology and trajectories and store it in ctx.obj
#Find more about angle calculations in functions.py

def main(ctx, top: str, trj: str):
    u = mda.Universe(top, trj)
    if 'ACE' in str(u.residues[0]):
        phi, psi = angle_to_list(u)
        ctx.obj['psi'] = psi
        ctx.obj['phi'] = phi
    else:
        phi, psi = prepare_dihedrals(top, u)
        ctx.obj['psi'] = psi
        ctx.obj['phi'] = phi

#Set up the commands
#Plot the data and give some informative visualizations

@main.command()
@click.pass_context
def plot(ctx):

    #Sets up the figure and subplots

    fig = plt.figure(constrained_layout=True)
    gs = gridspec.GridSpec(3, 2, figure=fig)

    #Scatter plot

    ax1 = fig.add_subplot(gs[0, 0])
    ax1.scatter(ctx.obj['phi'], ctx.obj['psi'], s=15, color='dimgray')
    ax1.plot([0, 0], [-180, 180], c='k', alpha=0.3)
    ax1.plot([-180, 180], [0, 0], c='k', alpha=0.3)
    ax1.set_xlim(-180, 180)
    ax1.set_ylim(-180, 180)
    ax1.set_xlabel('φ')
    ax1.set_ylabel('ψ')
    ax1.set_title('Ramachandran plot of a single amino acid')

    #Contour plot

    ax2 = fig.add_subplot(gs[0, 1])
    sns.kdeplot(ctx.obj['phi'], ctx.obj['psi'], ax=ax2, cmap='Reds', shade=True,shade_lowest=False)
    sns.rugplot(ctx.obj['phi'], color = 'k', ax=ax2)
    sns.rugplot(ctx.obj['psi'], color = 'k', vertical=True, ax=ax2)
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
    #Calculates angle differences

    psi_trans = angle_trans(ctx.obj['psi'])
    phi_trans = angle_trans(ctx.obj['phi'])
    psi_diff = angle_diff(psi_trans)
    phi_diff = angle_diff(phi_trans)

    combined_list = psi_diff + phi_diff

    #Angle difference over time for Phi

    ax3 = fig.add_subplot(gs[1, 0])
    #average_phi = sum(ctx.obj['phi'])/len(ctx.obj['phi'])
    #ax3.plot([0, len(ctx.obj['phi'])], [average_phi, average_phi], color='k', alpha=0.75, label='φ-avg: '+str(round(average_phi, 2)))
    #ax3.plot(range(0, len(ctx.obj['phi'])), ctx.obj['phi'], color='indianred')
    average_phi = sum(phi_diff) / len(phi_diff)
    stdev_phi = statistics.stdev(phi_diff)
    ax3.plot([0, len(phi_diff)], [average_phi, average_phi], color='k', alpha=0.75,
             label='φ-avg: '+str(round(average_phi, 2)) + ', stdev: ' + str(round(stdev_phi, 2)))
    ax3.plot(phi_diff, color='indianred', alpha=0.5)
    ax3.scatter(range(len(phi_diff)), phi_diff, marker = 'x', s = 5, color = 'k')
    ax3.set_ylabel('Angel difference')
    ax3.set_xlabel('Timeframe')
    ax3.set_title('φ-difference over time')
    ax3.set_ylim(min(combined_list)-5, max(combined_list)+5)
    ax3.set_xlim(0, len(ctx.obj['phi']))
    plt.legend(loc='upper right', borderaxespad=0.)

    #Angle difference over time for Psi

    ax4 = fig.add_subplot(gs[1, 1])
    #average_psi = sum(ctx.obj['psi'])/len(ctx.obj['psi'])
    #ax4.plot([0, len(ctx.obj['psi'])], [average_psi, average_psi], color='k', alpha=0.75, label='ψ-avg: '+str(round(average_psi, 2)))
    #ax4.plot(range(0, len(ctx.obj['psi'])), ctx.obj['psi'], color='darkkhaki')
    average_psi = sum(psi_diff) / len(psi_diff)
    stdev_psi = statistics.stdev(psi_diff)
    ax4.plot([0, len(psi_diff)], [average_psi, average_psi], color='k', alpha=0.75,
             label='ψ-avg: '+ str(round(average_psi, 2)) + ', stdev: ' + str(round(stdev_psi, 2)))
    ax4.plot(psi_diff, color='darkkhaki', alpha=0.5)
    ax4.scatter(range(len(psi_diff)), psi_diff, marker = 'x', s = 5, color = 'k')
    ax4.set_ylabel('Angel difference')
    ax4.set_xlabel('Timeframe')
    ax4.set_title('ψ-difference over time')
    ax4.set_ylim(min(combined_list)-5, max(combined_list)+5)
    ax4.set_xlim(0, len(ctx.obj['phi']))
    plt.legend(loc='upper right', borderaxespad=0.)

    bins_list_phi = list(range(int(min(ctx.obj['phi'])), int(max(ctx.obj['phi'])), 1))
    bins_list_psi = list(range(int(min(ctx.obj['psi'])), int(max(ctx.obj['psi'])), 1))

    #Histogram/Density plot for angle distribution

    ax5 = fig.add_subplot(gs[2,0:])
    #sns.kdeplot(phi_trans, ax=ax5, color = 'indianred', label='φ')
    #sns.kdeplot(psi_trans, ax=ax5, color = 'darkkhaki', label='ψ')
    #ax5.hist(ctx.obj['phi'], bins=bins_list_phi, alpha=0.75, color = 'indianred', label='φ')
    #ax5.hist(ctx.obj['psi'], bins=bins_list_psi, alpha=0.75, color = 'darkkhaki', label='ψ')
    sns.distplot(ctx.obj['phi'], ax=ax5, bins=bins_list_phi, color = 'indianred', label='φ')
    sns.distplot(ctx.obj['psi'], ax=ax5, bins=bins_list_psi, color = 'darkkhaki', label='ψ')
    ax5.set_xlim(-180, 180)
    ax5.set_title('Histogram/Density plot of angle distribution')
    plt.legend()

    plt.show()

#Print the data to stdout

@main.command()
@click.pass_context
def terminal(ctx):
    print(ctx.obj['psi'])
    print(ctx.obj['phi'])

def start():
    main(obj={})


if __name__ == '__main__':
    start()

