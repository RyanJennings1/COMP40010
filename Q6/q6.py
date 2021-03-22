#!/usr/bin/env python3
"""
  name: q6.py
  author: Ryan Jennings
  date: 2021-03-19
"""
from math import factorial
from typing import Dict, List, Union

import matplotlib.pyplot as plt
import numpy as np

from scipy.stats import poisson

LAMBDA: float = 0.003
MU: float = 1/433
S: int = 8
K: int = 10

RESULTS = Dict[str, Union[List[float], List[int], float]]

def main() -> None:
  """
  main method
  """
  #q6_1()
  #print_results(results=get_results_object())
  q6_4()

def q6_1() -> None:
  """
  Work out equations for Q6
  """
  n: List[int] = list(range(0, 26))
  for l in np.linspace(0.003, 0.018, num=16):
    Pns: List[float] = calc_Pns(n=n, _lambda=l)

    P0: float = 1/sum(Pns)
    Pn: List[float] = probabilities(n[1:], P0)
    Pn.insert(0, P0)
    Ns: List[float] = calc_Ns(n, Pn)
    N: float = sum(Ns)
    T: float = calc_T(N, P0)
    print(l, Pns[1], T)

def q6_4() -> None:
  """
  Apply modified Erlang B model to Q4
  """
  lambdas: List[float] = np.linspace(0.003, 0.018, num=16)
  Ws: List[float] = []
  Nqs: List[float] = []
  for l in lambdas:
    results = get_results_object(_lambda=l)
    Ws.append(results['W'])
    Nqs.append(results['Nq'])
  plot_Nqs_Ws(lambdas=lambdas, Nqs=Nqs, Ws=Ws)

def get_results_object(_lambda: float = LAMBDA) -> RESULTS:
  """
  Return an object containing all responses for a lambda
  """
  n: List[int] = list(range(0, 26))
  Pns: List[float] = calc_Pns(n, _lambda)
  P0: float = 1/sum(Pns)
  Pn: List[float] = probabilities(n[1:], P0, _lambda)
  Pn.insert(0, P0)
  Ns: List[float] = calc_Ns(n, Pn)
  N: float = sum(Ns)
  Nqs: List[float] = calc_Nqs(n, Pn)
  Nq: float = sum(Nqs)
  T: float = calc_T(N, P0, _lambda)
  W: float = calc_W(Nq, P0, _lambda)
  r: float = calc_r()
  return {
    'n': n, 'Pns': Pns, 'P0': P0, 'Pn': Pn, 'Ns': Ns,
    'N': N, 'Nqs': Nqs, 'Nq': Nq, 'T': T, 'W': W, 'r': r
  }

def print_results(results: RESULTS) -> None:
  """
  Print out all values
  """
  print(f"n: {results['n']}\n")
  print(f"Pns: {results['Pns']}\n")
  print(f"P0: {results['P0']}\n")
  print(f"Pn: {results['Pn']}\n")
  print(f"Ns: {results['Ns']}\n")
  print(f"N: {results['N']}\n")
  print(f"Nqs: {results['Nqs']}\n")
  print(f"Nq: {results['Nq']}\n")
  print(f"T: {results['T']}\n")
  print(f"W: {results['W']}\n")
  print(f"r: {results['r']}\n")

def probabilities(n: List[int], P0: float, _lambda: float = LAMBDA) -> List[float]:
  """
  Calculate probabilities for a range of n
  """
  arr: List[float] = []
  for i in n:
      if i > K:
          arr.append(0)
      else:
          if i >= min(S, K):
              x = ((_lambda/MU)**i) / (factorial(min(S, K)) * (min(S, K)**(i - min(S, K)))) * P0
              arr.append(x)
          else:
              #x = ((_lambda/MU)**i) / (factorial(i)*P0)
              x = (((_lambda/MU)**i) / (factorial(i)))*P0
              arr.append(x)
  return arr

def calc_Pns(n: List[int], _lambda: float = LAMBDA, K: float = K) -> List[float]:
  """
  Calculate probabilities for various n
  """
  arr: List[float] = []
  for i in n:
      if i <= min(S, K):
          x = ((_lambda/MU)**i) / factorial(i)
          arr.append(x)
      else:
          if i <= K:
              a = (_lambda/MU)**S
              b = (_lambda/(S*MU))**(i - S)
              x = (a * b) / factorial(S)
              arr.append(x)
          else:
              arr.append(0)
  return arr

def calc_Ns(n: List[int], Pn: List[float]) -> List[float]:
  """
  Calculate Number in system for various n and Pn
  """
  if len(n) != len(Pn):
      raise Exception("n and Pn are not the same length")
  return [n[c] * Pn[c] for c, i in enumerate(n)]

def calc_Nqs(n: List[int], Pn: List[float]) -> List[float]:
  """
  Calculate Number in queue for various n and Pn
  """
  if len(n) != len(Pn):
      raise Exception("n and Pn are not the same length")
  arr: List[float] = []
  for c, i in enumerate(n):
      if n[c] <= S:
          arr.append(0)
      else:
          arr.append((n[c] - S) * Pn[c])
  return arr

def calc_T(N: float, P0: float, _lambda: float = LAMBDA) -> float:
  """
  Calculate T
  """
  m = min(S, K)
  x = N / (_lambda*(1-(\
            ((_lambda/MU)**K) / (factorial(m)*(m**(K - m))) * P0)))
  return x

def calc_W(Nq: float, P0: float, _lambda: float = LAMBDA) -> float:
  """
  Calculate W
  """
  m = min(S, K)
  x = Nq / (_lambda*(1-(\
            ((_lambda/MU)**K) / (factorial(m)*(m**(K - m))) * P0)))
  return x

def calc_r(_lambda: float = LAMBDA) -> float:
  """
  Calculate r
  """
  return _lambda/(min(S, K)*MU)

def plot_Nqs_Ws(lambdas: List[float], Nqs: List[float], Ws: List[float]) -> None:
  """
  Plot the Number of requests in the queue and the
  queue wait time for K connections
  """
  fig, ax = plt.subplots(frameon=False)
  ax.set_axis_off()
  ax0 = fig.add_subplot(1, 2, 1)
  ax0.plot(lambdas, Nqs)
  ax0.set_title(f"Number in Queue vs Arrival Rate (s={S}, K={K})")
  ax0.set_ylabel('Nq (requests)')
  ax0.set_xlabel('λ (ms^-1)')

  ax1 = fig.add_subplot(1, 2, 2)
  ax1.plot(lambdas, Ws)
  ax1.set_title(f"Queue Wait Time vs Arrival Rate (s={S}, K={K})")
  ax1.set_ylabel('W (ms)')
  ax1.set_xlabel('λ (ms^-1)')
  plt.show()

#####################################################################
def scale(arr: List[float], a: int) -> List[float]:
  """
  Scale a non-numpy array
  """
  return list(map(lambda x: a*x, arr))

def erlang_c(arrival_rate: float, service_rate: float,
             K: int) -> float:
  """
  Erlang C formula
  K: int - number of connections
  Return probability of blocking
  """
  rho = arrival_rate * service_rate
  load = rho/K
  den: float = poisson.pmf(K, rho)
  return den / (den + (1-load)*(poisson.cdf(K-1, rho)))

def plot_Nq_W(K: int, lambdas: List[float]) -> None:
  """
  Plot the Number of requests in the queue and the
  queue wait time for K connections and a range of
  arrival rates
  """
  mu: int = MU
  Nqs: List[float] = []
  Ws: List[float] = []
  for l in lambdas:
    # Create data point
    rho = l * mu
    load = rho/K
    PD = erlang_c(l, mu, K)
    Nq = (PD*load) / (1-load)
    W = PD / ((1-load) * (K/mu))
    Nqs.append(Nq)
    Ws.append(W)
  fig, ax = plt.subplots(frameon=False)
  ax.set_axis_off()
  ax0 = fig.add_subplot(1, 2, 1)
  ax0.plot(lambdas, Nqs)
  ax0.set_title(f"Number in Queue vs Arrival Rate (K={K})")
  ax0.set_ylabel('Nq (requests)')
  ax0.set_xlabel('λ (ms^-1)')

  ax1 = fig.add_subplot(1, 2, 2)
  ax1.plot(lambdas, Ws)
  ax1.set_title(f"Queue Wait Time vs Arrival Rate (K={K})")
  ax1.set_ylabel('W (ms)')
  ax1.set_xlabel('λ (ms^-1)')
  plt.show()

if __name__ == "__main__":
    main()
