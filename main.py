from src.chemistry import Chemistry
from src.rocket import Rocket

chems = Chemistry.parse('test_real')

#print(chems)

rocket = Rocket(chems, 0.5, 1.1, 0.08)
