"""Every file filter of the desktop app's Open dialogs is one pywebview accepts (else its Choose button does nothing)."""
import re
from stagent.desktop import FILE_TYPES, FILE_FILTER

bad = [f for filters in FILE_TYPES.values() for f in filters if not re.search(FILE_FILTER, f)]
assert not bad, f"pywebview refuses these file filters: {bad}"
print(f"ok: {sum(map(len, FILE_TYPES.values()))} file filters")
