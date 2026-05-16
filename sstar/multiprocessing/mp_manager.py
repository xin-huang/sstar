# Copyright 2026 Xin Huang
#
# GNU General Public License v3.0
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, please see
#
#    https://www.gnu.org/licenses/gpl-3.0.en.html

import multiprocessing, queue, time
import numpy as np
from gaishi.generators import GenericGenerator
from multiprocessing import current_process, Manager, Process, Queue
from threading import Thread


def monitor(shared_dict: dict, workers: list[multiprocessing.Process]) -> None:
    """
    Monitors worker processes to ensure they complete successfully, initiating shutdown if any worker fails.

    Continuously checks the status of each worker process through a shared dictionary. If all workers
    have completed successfully, the monitoring loop exits. If any worker process terminates without
    marking its completion as 'Completed', a shutdown procedure is initiated for all workers.

    Parameters
    ----------
    shared_dict : dict
        A shared dictionary managed by a multiprocessing Manager. Worker processes use this dictionary
        to update their status. Keys are the names of worker processes, and values are status strings
        ('Completed', 'Failed', etc.).
    workers : list[multiprocessing.Process]
        A list of multiprocessing.Process objects, each representing a worker process to be monitored.

    Notes
    -----
    - The function assumes that worker processes update their status in the shared dictionary
      upon completion or failure.
    - In case of a worker failure (process is no longer alive but hasn't marked 'Completed'),
      `terminate_all_workers` is called to gracefully shutdown all workers.
    - The function uses a 1-second interval for periodic checks to balance responsiveness with efficiency.

    """
    while True:
        # alive_workers = [worker.name for worker in workers if worker.is_alive()]
        # completed_workers = [worker.name for worker in workers if shared_dict.get(worker.name) == 'Completed']
        # print("Monitoring", "Alive:", alive_workers, "Completed:", completed_workers)

        if all(shared_dict.get(worker.name) == "Completed" for worker in workers):
            # print("All workers completed their tasks successfully.")
            return

        for worker in workers:
            if not worker.is_alive() and shared_dict.get(worker.name) != "Completed":
                print(
                    f"{worker.name} did not complete successfully. Initiating shutdown."
                )
                terminate_all_workers(workers)
                print("All workers are terminated.")
                return

        time.sleep(1)  # Check periodically


def terminate_all_workers(workers: list[multiprocessing.Process]) -> None:
    """
    Terminates all worker processes and waits for them to complete.

    Sends a terminate signal to each worker process in the provided list and waits for each
    to join, ensuring all processes are properly terminated before proceeding.

    Parameters
    ----------
    workers : list[multiprocessing.Process]
        A list of multiprocessing.Process objects, each representing a worker process to be terminated.

    Notes
    -----
    - This function is typically called to ensure a clean shutdown in case of errors or when
      all work has been completed.
    - It first sends a `terminate` signal to each worker and then waits for each process to join,
      guaranteeing that no worker process is left hanging.

    """
    for w in workers:
        w.terminate()
    for w in workers:
        w.join()  # Wait for the process to terminate


def mp_manager(
    job: callable, data_generator: GenericGenerator, nprocess: int, **kwargs
) -> list:
    """
    Manages the distribution of tasks across multiple processes for parallel execution, ensuring
    reproducibility through controlled seed values for each task.

    This function initializes a pool of worker processes and distributes tasks among them.
    Each task involves executing a specified job, potentially with different seeds for
    each repetition to ensure variability yet reproducibility in stochastic processes.
    Results or errors are collected and returned as a list.

    Parameters
    ----------
    job : callable
        The function or callable object to be executed in parallel by each worker.
    data_generator : GenericGenerator
        An instance of a `GenericGenerator` subclass that yields dictionaries with parameters for each task.
        The `run` method in the corresponding job instance must be compatible with the parameters returned
        by the `get` method in the data_generator. This ensures that each task executed by the job function
        receives the correct parameters, facilitating reproducibility and consistency across tasks.
    nprocess : int
        The number of worker processes to use for executing the job in parallel. This determines
        the pool size of the multiprocessing environment.
    **kwargs : dict
        Additional keyword arguments to be passed directly to the job function. These are
        forwarded as-is to each job invocation.

    Returns
    -------
    res : list
        A list of collected results from each executed job.

    Raises
    ------
    Exception
        Captures and logs any exceptions encountered during the initialization or execution
        phase, including issues with starting worker processes or collecting results.

    Notes
    -----
    - The function utilizes a multiprocessing manager to create shared queues and dictionaries
      for task distribution and worker status monitoring.
    - To ensure smooth termination and cleanup, a monitoring thread is used to join all worker
      processes, and `cleanup_on_sigterm` is called to handle sudden terminations gracefully.

    """
    try:
        from pytest_cov.embed import cleanup_on_sigterm
    except ImportError:
        pass
    else:
        cleanup_on_sigterm()

    with Manager() as manager:
        res = []
        in_queue, out_queue = manager.Queue(), manager.Queue()
        shared_dict = manager.dict()
        workers = [
            Process(
                target=mp_worker, args=(in_queue, out_queue, shared_dict), kwargs=kwargs
            )
            for i in range(nprocess)
        ]

        for params in data_generator.get():
            in_queue.put((job, params))

        try:
            for w in workers:
                w.start()

            monitor_thread = Thread(target=monitor, args=(shared_dict, workers))
            monitor_thread.start()

            for i in range(len(data_generator)):
                items = out_queue.get()
                if items is None:
                    continue
                if isinstance(items, tuple) and "error" in items:
                    res = "error"
                    break
                res.extend(items)

            for w in workers:
                w.join()
        finally:
            for w in workers:
                w.terminate()
            monitor_thread.join()

    return res


def mp_worker(
    in_queue: queue.Queue, out_queue: queue.Queue, shared_dict: dict, **kwargs
) -> None:
    """
    A multiprocessing worker function that processes tasks from an input queue, executes a job,
    and reports the status to an output queue and a shared dictionary.

    This worker continuously fetches tasks from `in_queue`, each task comprising a repetition
    number, a seed value, and a job object. It executes the `run` method of the job object.
    Upon successful completion of a task, it places the result in `out_queue`. If the input queue
    is empty or an exception occurs during task processing, the worker updates its status in
    `shared_dict` and terminates gracefully.

    Parameters
    ----------
    in_queue : multiprocessing.managers.SyncManager.Queue
        The input queue from which the worker fetches tasks.
    out_queue : multiprocessing.managers.SyncManager.Queue
        The output queue where the worker posts the results of processed tasks.
    shared_dict : dict
        A shared dictionary where the worker updates its status. Uses the worker's process name
        as the key and the status ('Started', 'Completed', 'Failed') as the value.
    **kwargs : dict
        Additional keyword arguments that may be passed to the job's `run` method.

    Raises
    ------
    Exception
        Captures and logs any exceptions encountered during task processing. The worker updates
        its status as 'Failed' in `shared_dict` and posts an error message to `out_queue` before
        termination.

    Notes
    -----
    - The worker uses a timeout of 5 seconds for fetching tasks from `in_queue` to prevent
      indefinite blocking if the queue is empty.
    - Upon encountering an empty queue, the worker marks itself as 'Completed' in `shared_dict`
      and exits.
    - If an exception occurs, it marks itself as 'Failed' and posts the error to `out_queue`
      before breaking the loop and terminating.

    """
    process_name = current_process().name
    shared_dict[process_name] = "Started"

    while True:
        try:
            try:
                job, params = in_queue.get(timeout=5)
                items = job.run(**params, **kwargs)
            except queue.Empty:
                shared_dict[process_name] = "Completed"
                return  # Exit the loop and end the worker process
            out_queue.put(items)
        except Exception as e:
            shared_dict[process_name] = "Failed"
            out_queue.put(("error", str(e)))
            raise Exception(f"Worker {process_name} encountered an exception: {e}")
