"""Drop files from this node's page cache with posix_fadvise(DONTNEED).

Works without root, for clean pages only. Usage: python evict.py FILE...
"""

import os
import sys

os.sync()
for path in sys.argv[1:]:
    fd = os.open(path, os.O_RDONLY)
    try:
        os.posix_fadvise(fd, 0, 0, os.POSIX_FADV_DONTNEED)
    finally:
        os.close(fd)
