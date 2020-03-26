import MDAnalysis as mda
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from .functions import prepare_dihedrals, angle_to_list
import click

plt.style.use('seaborn-darkgrid')

@click.group()
@click.option('--topology', help='path to topology.gro')
@click.option('--xtc', help='path to trajectories.xtc')
@click.pass_context
def main(ctx, topology: str, xtc: str):
    u = mda.Universe(topology, xtc)
    if 'ACE' in str(u.residues[0]):
        phi, psi = angle_to_list(u)
        ctx.obj['psi'] = psi
        ctx.obj['phi'] = phi
    else:
        phi, psi = prepare_dihedrals(topology, u)
        ctx.obj['psi'] = psi
        ctx.obj['phi'] = phi


@main.command()
@click.pass_context
def plot(ctx):

    fig = plt.figure(constrained_layout=True)
    gs = gridspec.GridSpec(2, 2, figure=fig)
    ax1 = fig.add_subplot(gs[0, :])

    ax1.scatter(ctx.obj['phi'], ctx.obj['psi'], s=15, color='dimgray')
    ax1.plot([0, 0], [-180, 180], c='k', alpha=0.3)
    ax1.plot([-180, 180], [0, 0], c='k', alpha=0.3)
    ax1.set_xlim(-180, 180)
    ax1.set_ylim(-180, 180)
    ax1.set_xlabel('φ')
    ax1.set_ylabel('ψ')
    ax1.set_title('Ramachandran plot of a single amino acid')

    ax2 = fig.add_subplot(gs[1, 0])
    average_phi = sum(ctx.obj['phi'])/len(ctx.obj['phi'])
    ax2.plot([0, len(ctx.obj['phi'])], [average_phi, average_phi], color='k', alpha=0.75, label='φ-avg: '+str(round(average_phi, 2)))
    ax2.plot(range(0, len(ctx.obj['phi'])), ctx.obj['phi'], color='indianred')
    ax2.set_ylabel('Angel')
    ax2.set_xlabel('Timeframe')
    ax2.set_title('φ over time')
    ax2.set_ylim(-180, 180)
    ax2.set_xlim(0, len(ctx.obj['phi']))
    plt.legend()

    ax3 = fig.add_subplot(gs[1, 1])
    average_psi = sum(ctx.obj['psi'])/len(ctx.obj['psi'])
    ax3.plot([0, len(ctx.obj['psi'])], [average_psi, average_psi], color='k', alpha=0.75, label='ψ-avg: '+str(round(average_psi, 2)))
    ax3.plot(range(0, len(ctx.obj['psi'])), ctx.obj['psi'], color='darkkhaki')
    ax3.set_ylabel('Angel')
    ax3.set_xlabel('Timeframe')
    ax3.set_title('ψ over time')
    ax3.set_ylim(-180, 180)
    ax3.set_xlim(0, len(ctx.obj['phi']))
    plt.legend()

    plt.show()

@main.command()
@click.pass_context
def terminal(ctx):
    print(ctx.obj['psi'])
    print(ctx.obj['phi'])

def start():
    main(obj={})


if __name__ == '__main__':
    start()
