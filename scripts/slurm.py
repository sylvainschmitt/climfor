import dask
from dask_jobqueue.slurm import SLURMCluster

# dask
cluster = SLURMCluster(cores=10, memory='100GB')
cluster.adapt(minimum=1, maximum=10)  # Tells Dask to call `srun -n 1 ...` when it needs new workers

from dask.distributed import Client
client = Client(cluster)  # Connect this local process to remote workers

# init ee worker from https://gis.stackexchange.com/questions/480810/earth-engine-client-library-not-initialized-when-using-a-dask-cluster-and-xarray
class Plugin(dask.distributed.diagnostics.plugin.WorkerPlugin):
      def __init__(self, *args, **kwargs):
            pass  # the constructor is up to you
      def setup(self, worker: dask.distributed.Worker):
          ee.Initialize(opt_url='https://earthengine-highvolume.googleapis.com', project="silver-link-350720")
          pass
      def teardown(self, worker: dask.distributed.Worker):
          pass
      def transition(self, key: str, start: str, finish: str,
                      **kwargs):
          pass
      def release_key(self, key: str, state: str, cause: str | None, reason: None, report: bool):
          pass

plugin = Plugin()
client.register_plugin(plugin)
