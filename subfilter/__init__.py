
# Determine computation environment
try:
    print(f'Job id: {os.environ["SLURM_JOB_ID"]}')
    executing_on_cluster = True
except:
    print(f'Running interactively.')
    executing_on_cluster = False
