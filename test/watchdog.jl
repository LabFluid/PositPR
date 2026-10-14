# Execution guarded against infinite loops.
#
# Several PositPR functions use `while` loops whose stopping condition depends on
# the Posit arithmetic itself (eps, nextPosit, exp, log, sin, cos, pi, ...). If
# the arithmetic is wrong, the loop may never end and hang the whole suite.
#
# `guarded(f, args...)` runs `f(args...)` in a separate worker process
# (Distributed). If the time limit is exceeded, the process is killed, a new one
# is started on the next call and the case is reported as `TimedOut()`.

using Distributed

const WATCHDOG_TIMEOUT = 20.0          # seconds per call
const WATCHDOG_WORKER = Ref(0)

"Result of a call that did not finish within the time limit."
struct TimedOut end

function watchdog_worker()
    w = WATCHDOG_WORKER[]
    if w == 0 || !(w in workers())
        project = Base.active_project()
        w = only(addprocs(1; exeflags = `--project=$project`))
        Distributed.remotecall_eval(Main, w, :(using PositPR))
        WATCHDOG_WORKER[] = w
    end
    return w
end

function kill_worker(w)
    try
        kill(Distributed.worker_from_id(w).config.process)
    catch
    end
    try
        rmprocs(w; waitfor = 5)
    catch
    end
    WATCHDOG_WORKER[] = 0
end

"""
    guarded(f, args...; timeout = WATCHDOG_TIMEOUT)

Run `f(args...)` in a separate process. Returns the result, rethrows the
original exception (unwrapped from `RemoteException`) or returns `TimedOut()`.
"""
function guarded(f, args...; timeout = WATCHDOG_TIMEOUT)
    w = watchdog_worker()
    # Wait on a LOCAL task: querying the remote future of a hung worker would
    # hang as well.
    task = @async remotecall_fetch(f, w, args...)
    t0 = time()
    while !istaskdone(task)
        elapsed = time() - t0
        if elapsed > timeout
            kill_worker(w)
            return TimedOut()
        end
        elapsed < 0.05 ? yield() : sleep(0.005)
    end
    try
        return fetch(task)
    catch err
        while err isa TaskFailedException
            err = err.task.exception
        end
        err isa RemoteException && throw(err.captured.ex)
        throw(err)
    end
end

"Shut down the worker at the end of the suite."
stop_watchdog() = WATCHDOG_WORKER[] != 0 && kill_worker(WATCHDOG_WORKER[])
