"""공용 체로 가는 다리 — 경로 조작을 한 곳에 모은다.

`v3/code` 의 감사들은 `lib/` 을 import 하는 관례가 없다 (확인: `from
lib` 도 `import lib` 도 한 자리도 없다). 그런데 `mu` 를 인라인으로 짓는
스크립트가 303 개이고 (`v3/gate/SIEVE_INLINE`, `G95`), 그것들을 공용
몸통으로 옮기려면 각자 `sys.path` 를 만져야 한다 — 303 번 같은 조작을
복사하는 것은 **사본을 늘리지 말자**는 이 작업의 취지와 어긋난다.

그래서 조작은 여기 한 줄만 있고, 감사 쪽은 이렇게 부른다:

    from _sieve_shared import sieves        # (primes, Lambda, mu)
    from _sieve_shared import mu_upto       # mu 만 필요할 때

**이 파일은 구성을 하나도 담지 않는다.** 재수출뿐이라 `W8` 의 대조
집합이 안 늘어나고 (`def sieves*` 가 없다), `G95` 도 안 문다 (`mu` 에
부호를 뒤집는 자리가 없다). 몸통은 `lib/goldbach/sieve.py` 의 것 하나
이고 그것은 `audit_sieve.py` 의 것과 바이트 동일이며 `G16` 이 그
디렉터리도 훑는다.
"""
import os
import sys

# 몸통이 있는 뿌리를 **찾아서** 넣는다 — 칸을 세어 올라가지 않는다.
# `v3/code/` 에서는 두 칸 위가 뿌리지만 패킷에서는 `code/` 와 `lib/` 이
# 형제라 한 칸 위다. 칸 수를 박아 두면 패킷에서 `rounds/rNN` 을 뿌리로
# 잡고 **import 가 죽는다** — 실제로 죽어 있었고 r11 의 C 패스가 찾았다.
_HERE = os.path.dirname(os.path.abspath(__file__))
_REPO = None
_dir = _HERE
while True:
    if os.path.isfile(os.path.join(_dir, "lib", "goldbach", "sieve.py")):
        _REPO = _dir
        break
    _up = os.path.dirname(_dir)
    if _up == _dir:
        break
    _dir = _up
if _REPO is None:
    # 조용히 다른 몸통으로 넘어가지 않는다. 못 찾으면 선다.
    raise ImportError(
        "lib/goldbach/sieve.py 를 " + _HERE + " 위 어디에서도 못 찾았다")
if _REPO not in sys.path:
    sys.path.insert(0, _REPO)

from lib.goldbach.sieve import (  # noqa: E402
    mu_upto, primes_upto, sieves, spf_upto)

# **한계: `sieves(n)` 은 `n < 2^31` 에서만 옳다.** 몸통이 큰 소인수를
# 잡을 때 쓰는 `rem` 을 `np.arange(n + 1, dtype=np.int32)` 로 잡는데,
# `n + 1` 이 `int32` 를 넘으면 **numpy 가 예외를 안 내고 음수로 감는다**
# (`arange(2**31 + 11, dtype=int32)` 의 끝값이 `-2147483638` 이다).
# 그러면 `rem > 1` 이 거짓이 되어 그 자리의 `mu` 가 부호를 안 뒤집는다 --
# **조용히 틀린 `mu` 가 나온다.**
#
# `W8` 은 이것을 못 본다. 대조 범위가 `n <= 5e5` 라 한계 근처를 안 간다.
# **검사가 도는 범위와 코드가 쓰이는 범위가 다르다.**
#
# 이 한계를 몸통의 docstring 이 아니라 여기 적는 까닭: `_defbody` 가
# `def sieves(` 부터 들여쓰기 끝까지를 해시하므로 **docstring 한 줄이
# 해시를 바꾸고 `G16` 이 실패한다.** 정본을 손대려면 `SIEVE_HASHES` 도
# 같이 가야 하고 그건 별도 걸음이다.
#
# **실제 위반은 없다 (쟀다, 2026-09-04).** `v3/code` 에서 `2^31` 을 넘는
# 정수 리터럴을 전부 뽑아 용도를 봤다: `4294967296` 류는 모두 `say(...)`
# 안의 **산문**이고(「2^32 = 4294967296」을 인쇄한다), 유일한 상수
# `NCAP = 1000000000000` 은 `audit_dual_tailsweep.py` 의 **탐색 상한**
# (`while n <= NCAP`)이며 그 파일은 `mu` 체를 아예 안 돈다. 체 범위로
# 쓰이는 큰 수는 없다.
#
# **범위를 밝힌다:** 리터럴만 봤다. 범위를 **변수로 조립**하는 자리는
# 이 셈이 못 본다 -- 없는 것과 못 보는 것은 다르다. 그리고 배열 하나가
# `n` 바이트 이상을 먹으므로 `2^31` 근처는 물리적으로도 잘 안 간다는
# 것이 지금 위반이 없는 실질적인 까닭으로 보인다.

__all__ = ["sieves", "mu_upto", "primes_upto"]

# `primes_upto` 도 함께 내보내는 까닭: `W8` 은 몸통의 **본문만** 떼어
# 감사 스크립트 자신의 `primes_upto` 를 심고 돌린다. 그래서 몸통이
# 바이트 동일이어도 **다른 헬퍼 아래서는 다르게 행동하고 검사는 그
# 차이를 못 본다** (항목 48). `code/` 에 `primes_upto` 가 여섯 종
# 63 곳에 복제돼 있는데 `W8` 은 그중 어느 것도 안 본다.
# 그러므로 감사를 옮길 때 `mu` 체만 옮기고 헬퍼를 그 자리에 두면
# **옮기기 전보다 나빠질 수 있다** -- 몸통은 공용인데 헬퍼는 제각각인
# 상태가 된다. 여기서 둘을 함께 내보내 그 짝이 갈라지지 않게 한다.
