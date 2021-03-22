#!/usr/bin/env python3
"""
  name: q3.py
  author: Ryan Jennings
  date: 2021-03-13
"""
from typing import List

import matplotlib.pyplot as plt

LOAD: int = 0.5

def main() -> None:
  """
  main method
  """
  #plot_PL_vs_W()
  #plot_PL_vs_lambda()
  #plot_W_vs_lambda()
  q3_3()

def q3_3() -> None:
  """
  Plot values for Q3.3
  """
  pl: List[int] = [0.3333, 0.1429, 0.0667, 0.0323,
                   0.0159, 0.0079, 0.0039, 0.002]
  w: List[int] = [0, 16, 27, 34, 39, 43, 44, 46]
  fig = plt.figure()
  plt.title('Queue wait time vs Probability loss')
  plt.ylabel('W (ms)')
  plt.xlabel('PL (%)')
  plt.plot(pl, w)
  plt.show()

def plot_PL_vs_lambda() -> None:
  """
  Plot Probability Loss vs Arrival Rate
  """
  lambda_range: List[int] = [60, 80, 100, 120, 140, 160, 180]
  all_PLs: List[List[float]] = []
  for n in range(1, 9):
    PLs: List[float] = []
    for _lambda in lambda_range:
      l = 1/_lambda
      #m = 1/47
      m = 1/(2*_lambda)
      rho = load(l, m)
      PL = prob_loss(rho, n)
      PLs.append(PL)
    all_PLs.append(PLs)
    print(f"========================================")
    print(f"n = {n}")
    print(f"----------------------------------------")
    print(f"60\t80\t100\t120\t140\t160\t180\t")
    print(f"{PLs[0]:.4f}\t{PLs[1]:.4f}\t{PLs[2]:.4f}\t{PLs[3]:.4f}\t{PLs[4]:.4f}\t{PLs[5]:.4f}\t{PLs[6]:.4f}\t")
    print(f"========================================")
  fig, ax = plt.subplots(frameon=False)
  ax.set_axis_off()
  fig.tight_layout()
  for c, i in enumerate(all_PLs):
    ax0 = fig.add_subplot(2, 4, c+1)
    ax0.plot(lambda_range, all_PLs[c])
    ax0.set_title(f"PL vs 1/λ (n={c+1})")
    ax0.set_xlabel("1/λ")
    ax0.set_ylabel("PL")
    ax0.set_xlim([0, 180])
    ax0.set_ylim([0, 1])
  plt.show()

def plot_W_vs_lambda() -> None:
  """
  Plot Probability Loss vs Arrival Rate
  """
  lambda_range: List[int] = [60, 80, 100, 120, 140, 160, 180]
  all_Ws: List[List[float]] = []
  for n in range(1, 9):
    Ws: List[float] = []
    for _lambda in lambda_range:
      l = 1/_lambda
      #m = 1/47
      m = 1/(2*_lambda)
      rho = load(l, m)
      PL = prob_loss(rho, n)
      N = num_in_queue(rho, n)
      g = arrival_with_loss(l, PL)
      T = time_in_system(N, g)
      W = waiting_time(T, m)
      Ws.append(round(W))
    all_Ws.append(Ws)
    print(f"========================================")
    print(f"n = {n}")
    print(f"----------------------------------------")
    print(f"60\t80\t100\t120\t140\t160\t180\t")
    print(f"{Ws[0]}\t{Ws[1]}\t{Ws[2]}\t{Ws[3]}\t{Ws[4]}\t{Ws[5]}\t{Ws[6]}\t")
    print(f"========================================")
  fig, ax = plt.subplots(frameon=False)
  ax.set_axis_off()
  fig.tight_layout()
  for c, i in enumerate(all_Ws):
    ax0 = fig.add_subplot(2, 4, c+1)
    ax0.plot(lambda_range, all_Ws[c])
    ax0.set_title(f"W vs 1/λ (n={c+1})")
    ax0.set_xlabel("1/λ")
    ax0.set_ylabel("W")
    ax0.set_xlim([0, 180])
    ax0.set_ylim([0, 100])
  plt.show()

def plot_PL_vs_W() -> None:
  """Plot Probability Loss vs Queue Time"""
  lambda_range: List[int] = [60, 80, 100, 120, 140, 160, 180]
  all_PLs: List[List[float]] = []
  all_Ws: List[List[float]] = []
  for n in range(1, 9):
    PLs: List[float] = []
    Ws: List[float] = []
    for _lambda in lambda_range:
      l = 1/_lambda
      #m = 1/47
      m = 1/(2*_lambda)
      rho = load(l, m)
      PL = prob_loss(rho, n)
      N = num_in_queue(rho, n)
      g = arrival_with_loss(l, PL)
      T = time_in_system(N, g)
      W = waiting_time(T, m)
      PLs.append(PL)
      Ws.append(round(W))
    all_PLs.append(PLs)
    all_Ws.append(Ws)
    print(f"========================================")
    print(f"n = {n}")
    print(f"----------------------------------------")
    print(f"60\t80\t100\t120\t140\t160\t180\t")
    print(f"{PLs[0]:.4f}\t{PLs[1]:.4f}\t{PLs[2]:.4f}\t{PLs[3]:.4f}\t{PLs[4]:.4f}\t{PLs[5]:.4f}\t{PLs[6]:.4f}\t")
    print(f"{Ws[0]}\t{Ws[1]}\t{Ws[2]}\t{Ws[3]}\t{Ws[4]}\t{Ws[5]}\t{Ws[6]}\t")
    print(f"========================================")
  fig, ax = plt.subplots(frameon=False)
  ax.set_axis_off()
  fig.tight_layout()
  for c, i in enumerate(all_PLs):
    ax0 = fig.add_subplot(2, 4, c+1)
    ax0.plot(all_PLs[c], all_Ws[c])
    ax0.set_title(f"W vs PL (n={c+1})")
    ax0.set_xlabel("PL")
    ax0.set_ylabel("W")
    ax0.set_xlim([0, 1])
    ax0.set_ylim([0, 100])
  plt.show()

def load(_lambda: float, mu: float) -> float:
  """rho = lambda/mu"""
  return _lambda/mu

def prob_loss(load: float, n: int) -> float:
  """Probability of Loss (or Blocking) PL"""
  return ((1-load)*(load**n)) / (1-(load**(n+1)))

def num_in_queue(load: float, n: float) -> float:
  """Number of requests in the queue"""
  return (load+((load*n-n-1)*(load**(n+1)))) / ((1 - load**(n+1))*(1-load))

def arrival_with_loss(_lambda: float, pl: float) -> float:
  """Arrival rate with loss"""
  return _lambda * (1 - pl)

def time_in_system(N: float, gamma: float) -> float:
  """Time spent in system"""
  return N/gamma

def waiting_time(T: float, mu: float) -> float:
  """Waiting time for the queue"""
  return T - (1/mu)

if __name__ == "__main__":
    main()
