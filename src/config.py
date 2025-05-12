from pydantic import BaseModel
import numpy as np
import numpy.typing as npt
from enum import Enum
from toml import load
from pathlib import Path
from typing import Final


class Branch(Enum):
    Temporal = "temporal"
    Spatial = "spatial"


class ProblemType(Enum):
    Poiseuille = "Poiseuille"
    BoundaryLayer = "BoundaryLayer"
    Couette = "Couette"
    Custom = "Custom"


class ComplexNumber(BaseModel):
    r: float
    i: float

    def to_complex(self) -> complex:
        """Convert to a Python complex number."""
        return complex(self.r, self.i)


class VarsRange(BaseModel):
    min: float
    max: float
    num: int


class PlotLimits(BaseModel):
    xmin: float
    xmax: float
    ymin: float
    ymax: float

    # constructor based on the problem and branch
    @classmethod
    def default(cls, branch: Branch, problem: ProblemType) -> "PlotLimits":
        """Create default plot limits based on branch and problem type."""
        if branch == Branch.Temporal:
            if problem in {ProblemType.BoundaryLayer, ProblemType.Custom}:
                return cls(xmin=0.1, xmax=1.1, ymin=-1, ymax=0.1)
            else:
                return cls(xmin=0.1, xmax=1, ymin=-1, ymax=0.1)
        else:
            if problem in {ProblemType.BoundaryLayer, ProblemType.Custom}:
                return cls(xmin=-0.2, xmax=1, ymin=-0.1, ymax=0.8)
            else:
                return cls(xmin=0.2, xmax=1, ymin=-0.1, ymax=0.7)


class Config(BaseModel):
    DELTASTAR: Final[float] = 1.7207876573

    # params
    n: int
    re: float
    beta: complex

    # single run
    var: complex

    # multiple run
    vars_r: VarsRange
    vars_i: VarsRange
    vars_range_r: npt.NDArray[np.float64]
    vars_range_i: npt.NDArray[np.float64]

    # flags
    branch: Branch
    problem: ProblemType
    fileWriteEigenvalues: Path
    fileWriteEigenvector: Path
    doPlot: bool
    use_c: bool
    multipleRun: bool

    # custom problem flags
    filenameUprofile: Path
    plotUprofile: bool
    colX: int
    colY: int
    numSkipHeaderLines: int

    # plot
    plotLims: PlotLimits
    plotLabel: str

    class Config:
        arbitrary_types_allowed = True  # Allow arbitrary types like numpy.ndarray

    def __init__(self, **data):
        super().__init__(**data)  # Initialize normally

        # If problem is a BoundaryLayer, scale variables by the depth of delta*
        if self.problem == ProblemType.BoundaryLayer:
            scaling_factor = 1.0 / self.DELTASTAR
            self.re *= scaling_factor

            self.var *= scaling_factor
            self.beta *= scaling_factor

            self.vars_r.min *= scaling_factor
            self.vars_r.max *= scaling_factor

            self.vars_i.min *= scaling_factor
            self.vars_i.max *= scaling_factor

        self.vars_range_r = np.linspace(
            self.vars_r.min, self.vars_r.max, self.vars_r.num
        ).astype(np.float64)
        self.vars_range_i = np.linspace(
            self.vars_i.min, self.vars_i.max, self.vars_i.num
        ).astype(np.float64)

        if self.branch == Branch.Temporal:
            if self.use_c and np.abs(self.var) > 1e-10:
                self.plotLabel = "c"
            else:
                self.plotLabel = "omega"
        else:
            self.plotLabel = "alpha"

    @classmethod
    def from_toml(cls, file_path: str) -> "Config":
        """Loads configuration from a TOML file."""
        with open(file_path, "r") as file:
            data = load(file)

        # Convert complex numbers
        data["singleRunParams"]["var"] = ComplexNumber(**data["singleRunParams"]["var"]).to_complex()
        data["params"]["beta"] = ComplexNumber(**data["params"]["beta"]).to_complex()

        # flags
        data["flags"]["branch"] = Branch(data["flags"]["branch"])
        data["flags"]["problem"] = ProblemType(data["flags"]["problem"])
        data["flags"]["fileWriteEigenvalues"] = Path(
            data["flags"]["fileWriteEigenvalues"]
        )
        data["flags"]["fileWriteEigenvector"] = Path(
            data["flags"]["fileWriteEigenvector"]
        )

        # flags for custom problems
        data["customProblemFlags"]["filenameUprofile"] = Path(
            data["customProblemFlags"]["filenameUprofile"]
        )

        # default plot limits values
        default_plot_limits = PlotLimits.default(
            data["flags"]["branch"], data["flags"]["problem"]
        )

        # Flatten the structure for Pydantic
        parsed_data = {
            **data["params"],
            **data["singleRunParams"],
            **data["flags"],
            **data["customProblemFlags"],
            "vars_r": VarsRange(**data["multipleRunParams"]["vars_r"]),
            "vars_i": VarsRange(**data["multipleRunParams"]["vars_i"]),
            "plotLims": PlotLimits(**data["plot"]["plotLims"])
            if "plotLims" in data["plot"]
            else default_plot_limits,
            "vars_range_r": np.array([]),
            "vars_range_i": np.array([]),
            "plotLabel": "",
        }

        return cls(**parsed_data)  # Calls __init__, ensuring scaling if needed
