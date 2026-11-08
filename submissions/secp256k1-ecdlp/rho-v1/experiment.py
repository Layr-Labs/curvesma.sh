import json
import sys

def toy_rho(order):
    # Simple additive group simulation of rho collision
    # Returns True if collision found within expected steps
    seen = {}
    x = 1
    steps = 0
    max_steps = int(2 * (order ** 0.5))
    while steps < max_steps:
        if x in seen:
            return True
        seen[x] = steps
        x = (x * 2 + 1) % order  # simple map
        steps += 1
    return False

if __name__ == "__main__":
    data = json.load(sys.stdin)
    order = data["order"]
    success = toy_rho(order)
    print(json.dumps({"success": success, "order": order}))
