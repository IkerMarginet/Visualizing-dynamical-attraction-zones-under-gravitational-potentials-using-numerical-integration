
from __future__ import annotations

import numpy as np
import matplotlib.pyplot as plt
from dataclasses import dataclass
from typing import List, Tuple, Optional

# ---------------------------------------------------------------------------
# Gestion optionnelle de la barre de progression (tqdm)
# ---------------------------------------------------------------------------
try:
    from tqdm import tqdm
except ImportError:  # Fallback sans dépendance externe
    tqdm = None  # type: ignore


# ---------------------------------------------------------------------------
# Données : pôles / potentiels
# ---------------------------------------------------------------------------

@dataclass
class Attractor:
    """Représente un pôle attractif de type gravitationnel (V = -k/r)."""

    k: float  # Intensité k (>0) du potentiel  -k/r
    position: np.ndarray  # Coordonnées (x, y) fixes du pôle
    color: Tuple[float, float, float]  # Couleur RGB ∈ [0, 1]³ utilisée pour la carte


# ---------------------------------------------------------------------------
# Force résultante exercée sur la particule
# ---------------------------------------------------------------------------

def force_on_particle(pos: np.ndarray, attractors: List[Attractor]) -> np.ndarray:
    """Calcule la force totale **F**(pos) = Σ –k r̂ / r² exercée par les pôles."""

    f = np.zeros(2)
    for att in attractors:
        r_vec = pos - att.position
        r = np.linalg.norm(r_vec)
        if r < 1e-9:  # Singularity guard when the particle meets the pole
            continue
        f -= att.k * r_vec / r**3
    return f


# ---------------------------------------------------------------------------
# Intégrateurs temporels
# ---------------------------------------------------------------------------

def integrate_rk4(
    pos0: np.ndarray,
    vel0: np.ndarray,
    attractors: List[Attractor],
    dt: float,
    n_steps: int,
    r_stop: float,
) -> Optional[int]:
    """Intégrateur Runge–Kutta 4ᵉ ordre.

    Retourne l’indice du pôle capturant ou *None* si la particule s’échappe.
    """

    pos = pos0.copy()
    vel = vel0.copy()

    def acc(p):
        return force_on_particle(p, attractors)

    for _ in range(n_steps):
        k1v = acc(pos) * dt
        k1p = vel * dt

        k2v = acc(pos + 0.5 * k1p) * dt
        k2p = (vel + 0.5 * k1v) * dt

        k3v = acc(pos + 0.5 * k2p) * dt
        k3p = (vel + 0.5 * k2v) * dt

        k4v = acc(pos + k3p) * dt
        k4p = (vel + k3v) * dt

        vel += (k1v + 2 * k2v + 2 * k3v + k4v) / 6.0
        pos += (k1p + 2 * k2p + 2 * k3p + k4p) / 6.0

        # Conditions d’arrêt
        for idx, att in enumerate(attractors):
            if np.linalg.norm(pos - att.position) < r_stop:
                return idx  # Capturé
        if np.linalg.norm(pos) > 2.0:
            return None  # Échappé

    return None  # Temps max atteint sans capture


def integrate_symplectic(
    pos0: np.ndarray,
    vel0: np.ndarray,
    attractors: List[Attractor],
    dt: float,
    n_steps: int,
    r_stop: float,
) -> Optional[int]:
    """Intégrateur d’Euler symplectique (ordre 1, conservatif)."""

    pos = pos0.copy()
    vel = vel0.copy()

    for _ in range(n_steps):
        # Demi‑pas vitesse
        vel += 0.5 * force_on_particle(pos, attractors) * dt
        # Pas complet position
        pos += vel * dt
        # Demi‑pas vitesse restant
        vel += 0.5 * force_on_particle(pos, attractors) * dt

        # Conditions d’arrêt
        for idx, att in enumerate(attractors):
            if np.linalg.norm(pos - att.position) < r_stop:
                return idx
        if np.linalg.norm(pos) > 2.0:
            return None

    return None


# ---------------------------------------------------------------------------
# Génération de la carte (progress bar incluse)
# ---------------------------------------------------------------------------

def generate_map(
    attractors: List[Attractor],
    *,
    integrator: str = "rk4",
    grid_size: int = 400,
    dt: float = 0.01,
    n_steps: int = 5000,
    r_stop: float = 0.05,
    v0: np.ndarray | None = None,
    show_progress: bool = True,
) -> np.ndarray:
    """Calcule une *image* RGB (H×W×3) représentant la cartographie.

    Le paramètre `show_progress` active (par défaut) une barre de progression
    *tqdm* si la librairie est installée.
    """

    if v0 is None:
        v0 = np.zeros(2)

    # Choix de l’intégrateur
    if integrator == "rk4":
        integrator_func = integrate_rk4
    elif integrator in {"symplectic", "euler_symplectic"}:
        integrator_func = integrate_symplectic
    else:
        raise ValueError(f"Intégrateur inconnu : {integrator}")

    result_indices = np.full((grid_size, grid_size), -1, dtype=int)

    xs = np.linspace(-1.0, 1.0, grid_size)
    ys = np.linspace(-1.0, 1.0, grid_size)

    total_px = grid_size * grid_size
    iterator = range(total_px)

    # Activation éventuelle de la barre de progression
    if show_progress and tqdm is not None:
        iterator = tqdm(iterator, desc="Calcul de la carte", unit="px")

    # Boucle unique pour bénéficier d’une barre % correcte
    for p in iterator:
        i = p // grid_size  # ligne (y)
        j = p % grid_size   # colonne (x)
        x = xs[j]
        y = ys[i]

        idx = integrator_func(np.array([x, y]), v0, attractors, dt, n_steps, r_stop)
        result_indices[i, j] = -1 if idx is None else idx

    # Traduction indices → couleurs
    img = np.zeros((grid_size, grid_size, 3))
    white = (1.0, 1.0, 1.0)
    for idx, att in enumerate(attractors):
        img[result_indices == idx] = att.color
    img[result_indices == -1] = white  # Trajectoires échappées

    return img


# ---------------------------------------------------------------------------
# Exemple complet avec barre de progression
# ---------------------------------------------------------------------------

def demo():  # noqa: D401 – simple demo
    """Génère et affiche la cartographie pour trois pôles avec deux intégrateurs."""

    attractors = [
        Attractor(k=1.0, position=np.array([-0.5, 0.0]), color=(1.0, 0.0, 0.0)),
        Attractor(k=1.0, position=np.array([0.5, 0.0]), color=(0.0, 1.0, 0.0)),
        Attractor(k=1.0, position=np.array([0.0, 0.8]), color=(0.0, 0.0, 1.0)),
    ]

    # Runge–Kutta 4
    img_rk4 = generate_map(
        attractors,
        integrator="rk4",
        grid_size=450,
        dt=0.005,
        n_steps=4000,
        show_progress=True,
    )

    # Euler symplectique
    img_symp = generate_map(
        attractors,
        integrator="symplectic",
        grid_size=450,
        dt=0.005,
        n_steps=4000,
        show_progress=True,
    )

    # Affichage côte‑à‑côte
    fig, axes = plt.subplots(1, 2, figsize=(12, 6), dpi=100)
    axes[0].imshow(img_rk4, extent=[-1, 1, -1, 1], origin="lower")
    axes[0].set_title("Runge–Kutta 4ᵉ ordre")

    axes[1].imshow(img_symp, extent=[-1, 1, -1, 1], origin="lower")
    axes[1].set_title("Euler symplectique")

    for ax in axes:
        ax.set_xlabel("x")
        ax.set_ylabel("y")
        ax.set_aspect("equal", adjustable="box")

    plt.tight_layout()
    plt.show()



