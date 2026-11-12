"""Local dask cluster for pbzarr computations, after pbview's cluster.py."""

from contextlib import contextmanager

from dask.distributed import Client, LocalCluster


class DaskClient:
    def __init__(self, workers=2, threads=2, memory_limit="auto"):
        self.n_workers = workers
        self.threads_per_worker = threads
        self.memory_limit = memory_limit

    @contextmanager
    def get_client(self):
        cluster = LocalCluster(
            n_workers=self.n_workers,
            threads_per_worker=self.threads_per_worker,
            memory_limit=self.memory_limit,
        )
        client = Client(cluster)
        try:
            yield client
        finally:
            client.close()
            cluster.close()
