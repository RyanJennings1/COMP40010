#!/usr/bin/env python3
"""
  name: q4.py
  author: Ryan Jennings
  date: 2021-03-15
"""
from typing import List

import matplotlib.pyplot as plt
import numpy as np

from scipy.stats import poisson

MU: int = 433

def main() -> None:
  """
  main method
  """
  q4_4()

def q4_1() -> None:
  """
  Plot values for Q4.1
  """
  _lambda: List[float] = np.linspace(0.003, 0.023, num=21)
  plot_Nq_W(K=10, lambdas=_lambda)

def q4_2() -> None:
  """
  Plot values for Q4.2
  """
  _lambda: List[float] = np.linspace(0.003, 0.018, num=16)
  plot_Nq_W(K=8, lambdas=_lambda)

def q4_3() -> None:
  """
  Plot values for Q4.3

  Subplots:
  k=6  k=4  k=2
  L_Nq L_Nq L_Nq
  L_W  L_W  L_W
  """
  _lambda: List[float] = np.linspace(0.003, 0.018, num=16)
  #mu: int = 26
  mu: int = MU
  all_Nqs: List[List[float]] = []
  all_Ws: List[List[float]] = []
  Ks: List[int] = [6, 4, 2]
  for K in Ks:
    Nqs: List[float] = []
    Ws: List[float] = []
    for l in _lambda:
      # Create data point
      rho = l * mu
      load = rho/K
      PD = erlang_c(l, mu, K)
      Nq = (PD*load) / (1-load)
      W = PD / ((1-load) * (K/mu))
      Nqs.append(Nq)
      Ws.append(W)
    all_Nqs.append(Nqs)
    all_Ws.append(Ws)
  fig, ax = plt.subplots(frameon=False)
  ax.set_axis_off()
  #fig.tight_layout()
  for c, i in enumerate(all_Nqs):
    ax0 = fig.add_subplot(2, 3, c+1)
    ax0.plot(_lambda, all_Nqs[c])
    ax0.set_title(f"Number in Queue vs Arrival Rate, K={Ks[c]}")
    ax0.set_ylabel('Nq (requests)')
    ax0.set_xlabel('λ (ms^-1)')
  for c, i in enumerate(all_Ws):
    ax1 = fig.add_subplot(2, 3, c+1+3)
    ax1.plot(_lambda, all_Ws[c])
    ax1.set_title(f"Queue Wait Time vs Arrival Rate, K={Ks[c]}")
    ax1.set_ylabel('W (ms)')
    ax1.set_xlabel('λ (ms^-1)')
  plt.show()

def scale(arr: List[float], a: int) -> List[float]:
  """
  Scale a non-numpy array
  """
  return list(map(lambda x: a*x, arr))

def q4_4() -> None:
  """
  Plot values for Q4.4
  """
  #_lambda: List[float] = np.linspace(0.003, 0.018, num=16)
  a: int = 4
  _lambda: List[float] = scale([0.003, 0.006, 0.012, 0.024, 0.048], 2)
  #mu: int = 26
  mu: int = MU
  Ks: List[int] = scale([1, 2, 4, 8, 16], a)
  Nqs: List[float] = []
  Ws: List[float] = []
  for i, l in enumerate(_lambda):
    # Create data point
    K = Ks[i]
    rho = l * mu
    load = rho/K
    PD = erlang_c(l, mu, K)
    Nq = (PD*load) / (1-load)
    W = PD / ((1-load) * (K/mu))
    Nqs.append(Nq)
    Ws.append(W)
  fig, ax = plt.subplots(frameon=False)
  ax.set_axis_off()
  #fig.tight_layout()
  ax0 = fig.add_subplot(1, 2, 1)
  ax0.plot(_lambda, Nqs)
  ax0.set_title(f"Number in Queue vs Arrival Rate (Doubling)")
  ax0.set_ylabel('Nq (requests)')
  ax0.set_xlabel('λ (ms^-1)')

  ax1 = fig.add_subplot(1, 2, 2)
  ax1.plot(_lambda, Ws)
  ax1.set_title(f"Queue Wait Time vs Arrival Rate (Doubling)")
  ax1.set_ylabel('W (ms)')
  ax1.set_xlabel('λ (ms^-1)')
  plt.show()

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
