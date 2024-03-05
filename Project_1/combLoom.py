import loompy

loompy.combine(['/home/flyhunter/Wang/data/control/velocyto/control.loom', '/home/flyhunter/Wang/data/stress/velocyto/stress.loom'], key="Accession", output_file='controlStress.loom')
