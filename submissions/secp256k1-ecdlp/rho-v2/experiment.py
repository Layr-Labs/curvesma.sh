import json
import sys
import random

def toy_rho(order, seed=1):
    random.seed(seed)
    # Additive cyclic group simulation
    # Recover discrete log of target by rho-style walk
    target = random.randrange(order)
    # Simple map: x -> 2x + 1 mod order (for illustration)
    def step(x):
        return (2 * x + 1) % order
    seen = {}
    x = random.randrange(order)
    steps = 0
    max_steps = int(4 * (order ** 0.5)) + 10
    while steps < max_steps:
        if x in seen:
            # Collision found; in a real ECDLP we would solve for the log
            return {"success": True, "steps": steps, "order": order}
        seen[x] = steps
        x = step(x)
        steps += 1
    return {"success": False, "steps": steps, "order": order}

if __name__ == "__main__":
    data = json.load(sys.stdin)
    order = data["order"]
    result = toy_rho(order)
    print(json.dumps(result))
