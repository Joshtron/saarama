import MDAnalysis as mda
import matplotlib.pyplot as plt
from .functions import calc_q_vectors, calc_cross_vectors, \
    calc_normals, calc_orthogonal_unit_vectors, calc_dihedral_angle, \
    get_dihedrals, prepare_dihedrals
import click

@click.group()
@click.option('--topology', help='path to topology.gro')
@click.option('--xtc', help='path to trajectories.xtc')
@click.pass_context
def main(ctx, topology: str, xtc: str):
    u = mda.Universe(topology, xtc)
    psi, phi = prepare_dihedrals(topology, u)
    ctx.obj['psi'] = psi
    ctx.obj['phi'] = phi

@main.command()
@click.pass_context
def plot(ctx):
    plt.scatter(ctx.obj['psi'], ctx.obj['phi'], s=15)
    plt.xlim(-180, 180)
    plt.ylim(-180, 180)
    plt.xlabel('φ')
    plt.ylabel('ψ')
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
