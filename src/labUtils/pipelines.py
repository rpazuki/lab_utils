##############################################################
#  labUtils: A python package for laboratory data analysis   #
#  and automation.                                           #
#                                                            #
#  Author: Roozbeh H. Pazuki                                 #
#  License: MIT                                              #
##############################################################
import importlib
import logging
from abc import ABC, abstractmethod
from pathlib import Path
from typing import Any, Callable, Mapping, Optional, TypeVar

import yaml
from addict import Dict as DefaultDict
from pandas import DataFrame

_Self = TypeVar('_Self', bound='AbstractPipeline')


# logger
log = logging.getLogger(__name__)


class Dict(DefaultDict):
    def __missing__(self, key) -> None:
        # raise KeyError(key)
        # calling dict.unassinged properties return None
        return None

class IncompatibleArgsException(Exception):
    """Exception raised when the arguments passed to a process are incompatible.
    """
    def __init__(self, message):
        self.message = message
        super().__init__(self.message)


class AbstractPipeline(ABC):
    """An abstract pipeline class.
    """

    def __init__(self, processes=[]):
        self.processes = processes

    def append_process(self, process) -> None:
        """Append a process to the pipeline."""
        self.processes.append(process)

    def __call__(self, /, **kwargs) -> Dict:
        return self.process(**kwargs)

    @abstractmethod
    def process(self, /, **kwargs) -> Dict:
        """Process and return the pipeline."""

    @abstractmethod
    def __rshift__(self: _Self, other) -> _Self:
        pass

    @abstractmethod
    def __mul__(self, other) -> 'ProcessFork':
        pass


class AbstractProcess(ABC):

    @abstractmethod
    def __call__(self, **kwargs) -> Dict:
        """The operation of the process must happen here.

        It can define any number of arguments it like and
        must have one '**kwargs' at the end.

        The returns will be named payload as a Dict object
        and will be passed to the down-stream process or
        return to the caller.
        """


    def __rshift__(self, other) -> AbstractPipeline:
        """Appends the process to the end of an Process or DFPipeline.

           This is an imumutable operation and new DFPipeline will be
           returned.

           Example:
                proc1 >> proc2 >> proc3

        Parameters
        ----------
        other : AbstractProcess | DFPipeline
            The process or pipeline that will be appended to the given
            process.
            For ProcessFork in the RHS, the output of the ProcessFork
            will be joined by returning a ProcessJoined object.

        Returns
        -------
        AbstractPipeline
            A new (immutable) pipeline object that contains the current
            process followed by 'other' argument.

        Raises
        ------
        ValueError
            Arises when the RHS is not one of the accptable types.
        """
        if issubclass(type(other), AbstractProcess):
            if isinstance(other, ProcessFork):
                other = ProcessJoined(other)
            return DFPipeline([self, other])
        elif issubclass(type(other), AbstractPipeline):
            return DFPipeline(
                [self] + other.processes)
        else:
            raise ValueError(
                f"The '{type(other)}' must be a AbstractProcess or DFPipeline.")

    def __mul__(self, other):
        """Fork two or more processes.

           This is an imumutable operation and new ProcessFork will be
           returned.

           Example:
               p1 * p2 : Forks two processes
               p1 * p2 * p3 : Forks three processes

        Parameters
        ----------
        other : int, ProcessFork, AbstractProcess, DFPipeline
            When the RHS is n (int), it is forked to n replicate.
            When the RHS is ProcessFork, AbstractProcess or DFPipeline,
            the LHS process stack on its process.

        Returns
        -------
        ProcessFork
            A callable object that forks the processes.

        Raises
        ------
        ValueError
            Raises when the 'other' argument is not an int, a ProcessFork, sub class of
            AbstractProcess, or DFPipeline.
            Also, in case for RHS int and LHS ProcessFork, it will raise the error.
        """
        if isinstance(other, int):
            if isinstance(self, ProcessFork):
                raise ValueError(
                    "The 'ProcessFork' has already forked (cannot be multiplied).")
            else:
                return ProcessFork([self] * other)
        if isinstance(other, ProcessFork):
            if isinstance(self, ProcessFork):
                return ProcessFork(self.processes + other.processes)
            else:
                return ProcessFork([self] + other.processes)
        if issubclass(type(other), AbstractProcess):
            return ProcessFork([self, other])
        if issubclass(type(other), DFPipeline):
            return ProcessFork([self, other])
        else:
            raise ValueError(
                f"The '{type(other)}' must be an int, a AbstractProcess, DFPipeline "
                f"or tuples of both.")


class Process(AbstractProcess):
    """A singleton object for instantiable processes.
    """
    _instance = None

    def __new__(cls):
        if cls._instance is None:
            cls._instance = super(Process, cls).__new__(cls)
            # Put any initialization here.
        return cls._instance


class ProcessPassThrough(Process):
    def __call__(self, **kwargs) -> Dict:
        """Payload is simplely pass to the next process."""
        return Dict(**kwargs)


class ProcessJoined(AbstractProcess):
    def __init__(self,
                 forkedProcess: AbstractProcess,
                 kwargs_mapping: Mapping[int, list[tuple[str, str]]] = {}):
        """Join all the forked processes and returns as a single Mapping.

        Parameters
        ----------
        forkedProcess : Process
            The process that has been forked.
        kwargs_mapping : Mapping[int, list[tuple[str, str]]], optional
            Provides the renaming of the processes outputs, by default None.
            The key is the index of the process in the forked process.
            The value is a list of tuples of the old key and the new key.
        """
        self.forkedProcess = forkedProcess
        self.kwargs_mapping = kwargs_mapping

    def __call__(self, **kwargs) -> Dict:
        # First, calls all the sub-processes in the fork.
        tuple_return = self.forkedProcess(**kwargs)
        # Next, if there is any renaming requested in the
        # payloads of the forked sub-processes, we do it here
        for index, names in self.kwargs_mapping.items():
            ret = tuple_return[index]
            for old_key, new_key in names:
                ret[new_key] = ret.pop(old_key)
        # Finally, make a union of all the retuned payloads.
        # Note: the names from returns of higher rank  has precedence
        new_kwargs = tuple_return[0]
        for d in tuple_return[1:]:
            new_kwargs |= d
        return new_kwargs


class ProcessFork(AbstractProcess):
    def __init__(self,
                 processes: list):
        assert len(processes) > 1, "The processes must be more than one."
        self.processes = processes

    def __getitem__(self, key: int) -> AbstractProcess | AbstractPipeline:
        return self.processes[key]

    def __setitem__(self, key: int, process) -> None:
        self.processes[key] = process

    def __call__(self,  **kwargs) -> tuple[Dict, ...]:
        """It calles each sub-processes of the fork and returns thier payload as a tuple."""
        return tuple(p(**kwargs) for p in self.processes)

    def __rshift__(self, other) -> ProcessJoined:
        # It first join the forked instance, and next,
        # calls the abstractProcess rshift to append the
        # RHS to the joined process.
        return ProcessJoined(self) >> other

    def __itruediv__(self, other) -> ProcessJoined:
        return self.__truediv__(other)

    def __truediv__(self, other: Mapping[int, list[tuple[str, str]]]) -> ProcessJoined:
        """combine the output of a forked process.

            Examples:
            forked_process /= {0:[ ('old_name1', 'new_name1'),
                                   ('old_name2', 'new_name2')],
                               1:[ ('old_name3', 'new_name3'),
                                   ('old_name4', 'new_name4')]
                                   }

        Parameters
        ----------
        other : Mapping[int, list[tuple[str, str]]]
            A mapping of .

        Returns
        -------
        ProcessJoined
            A joined processes.

        Raises
        ------
        ValueError
            raises when the 'other' argument is not a type of a mapping.
        """
        if isinstance(other, Mapping):
            error_msg = ("The number of renaming must be less than "
                      "or equal to the number of processes.")
            assert len(other) <= len(self.processes), error_msg
            return ProcessJoined(self, kwargs_mapping=other)
        else:
            raise ValueError(
                f"The {type(other)=} must be a Mapping[str, int] ")


class DFPipeline(AbstractPipeline):

    def __rshift__(self, other) -> AbstractPipeline:
        """Appends the process to the end of an Process or DFPipeline.

           This is an imumutable operation and new DFPipeline will be
           returned.

           Example:
                proc1 >> proc2 >> proc3

        Parameters
        ----------
        other : AbstractProcess | DFPipeline
            The process or pipeline that will be appended to the given
            pipeline.
            For ProcessFork in the RHS, the output of the ProcessFork
            will be joined by returning a ProcessJoined object inside a
            pipeline.

        Returns
        -------
        AbstractPipeline
            A new (immutable) pipeline object that contains the current
            pipeline followed by 'other' argument.

        Raises
        ------
        ValueError
            Arises when the RHS is not one of the accptable types.
        """
        if issubclass(type(other), AbstractProcess):
            if isinstance(other, ProcessFork):
                other = ProcessJoined(other)
            return DFPipeline(self.processes + [other])
        elif isinstance(other, DFPipeline):
            return DFPipeline(self.processes + other.processes)
        else:
            raise ValueError(
                f"The '{type(other)}' must be a Process or DFPipeline.")

    def __mul__(self, other) -> ProcessFork:
        """Fork two or more pipslines.

           This is an imumutable operation and new ProcessFork will be
           returned.

           Example:
               p1 * p2 : Forks two processes/pipelines
               p1 * p2 * p3 : Forks three processes/piplines

        Parameters
        ----------
        other : int, ProcessFork, AbstractProcess, DFPipeline
            When the RHS is n (int), it is forked to n replicate.
            When the RHS is ProcessFork, AbstractProcess or DFPipeline,
            the LHS process stack on its process.

        Returns
        -------
        ProcessFork
            A callable object that forks the processes.

        Raises
        ------
        ValueError
            Raises when the 'other' argument is not an int, a ProcessFork, sub class of
            AbstractProcess, or DFPipeline.
            Also, in case for RHS int and LHS ProcessFork, it will raise the error.
        """
        if isinstance(other, int):
            return ProcessFork([self] * other)
        if isinstance(other, ProcessFork):
            raise ValueError(
                f"The {type(other)=} cannot be mutiplied from LHS of a DFPipeline."
            )
        if issubclass(type(other), AbstractProcess):
            return ProcessFork([self, other])
        if issubclass(type(other), DFPipeline):
            return ProcessFork([self, other])
        else:
            raise ValueError(
                f"The {type(other)=} must be an int, a Process, DFPipeline "
                f"or tuples of both.")

    def process(self, /, **kwargs) -> Dict:
        """Proccess and return the pipeline.

        all kwargs are default values for any processes in down-stream. However,
        the process returned parameters have precedence.

        Returns
        -------
        Dict
            A dict object that conains the processed dtataframe and other parameters.

        Raises
        ------
        IncompatibleArgsException
        """
        index = 0
        try:
            payload_kwargs = Dict(**kwargs)
            for index, process in enumerate(self.processes):
                ret: Dict = process(**payload_kwargs)
                # Union the returned payload with previous ones.
                # This will be passed to next process or return to the caller.
                #
                #  1- Process parameters have precedence over the kwargs.
                #  2- The latest process parameters have precedence over the formeres.
                payload_kwargs |= ret
        except TypeError as e:
            # The name of the class
            def name(obj) -> str:
                return type(obj).__name__

            if (len(e.args) > 0 and
                    ("missing 1 required keyword-only argument" in e.args[0] or
                     "missing 1 required positional argument" in e.args[0]
                     )):
                if index > 0:
                    previous_process = name(self.processes[index - 1])
                else:
                    previous_process = "(input of the pipline)"
                raise IncompatibleArgsException(
                    f"The process '{name(process)}' received incompatible payload from "
                    f"the previous process '{previous_process}'") from e
            else:
                raise e
        return payload_kwargs


class ProcessLogic(AbstractProcess):
    def __init__(self, logic_callback: Callable[..., Dict]):
        self.logic_callback = logic_callback

    def __call__(self, **kwargs) -> Dict:
        return self.logic_callback(**kwargs)


class ProcessLogicProperty(AbstractProcess):
    def __init__(self, logic_callback: Callable[..., Dict]):
        self.logic_callback = logic_callback
        self.caller_class = None

    def __call__(self, **kwargs) -> Dict:
        return self.logic_callback(self.caller_class, **kwargs)


class ProcessFactory:
    def __init__(self, factory_callback: Callable[..., AbstractProcess]):
        self.factory_callback = factory_callback

    def __call__(self, *args: Any, **kwargs) -> AbstractProcess:
        return self.create(*args, **kwargs)

    def create(self, *args: Any, **kwargs) -> AbstractProcess:
        """Create a parametrised process (arguments) to be attached and called later."""
        process = self.factory_callback(*args, **kwargs)
        return process


class InputProcess(Process):
    def __call__(self, *, name:str, src:str, package:str, method:str, **kwargs) -> Dict:
        """Input process to load data from different sources.

        Parameters
        ----------
        name : str
            The name of the data source.
        src : str
            The source path or URL of the data.
        package : str
            The package name to use for loading the data.
        method : str
            The method name within the package to use for loading the data.
        """
        pkg = importlib.import_module(package)
        func = getattr(pkg, method)
        data = func(Path(src), **kwargs)
        return Dict({f"{name}":data})


class DFProcess(Process):
    def __call__(self, *, payload:Dict, name:str, package:str, method:str, parameters:dict, **kwargs) -> Dict:
        """DataFrame process to save data to different destinations.

        Parameters
        ----------
        payload : Dict
            The data payload to be saved.
        name : str
            The name of the DataFrame.
        package : str
            The package name to use for saving the DataFrame.
        method : str
            The method name within the package to use for saving the DataFrame.
        parameters:dict

        """
        def is_hashable(obj):
            try:
                hash(obj)
                return True
            except TypeError:
                return False
        pkg = importlib.import_module(package)
        func = getattr(pkg, method)
        arguments = {}
        for key, value in parameters.items():
            if is_hashable(value) and value in payload:
                arguments[key] = payload[value]
            else:
                arguments[key] = value
        data = func(**arguments)
        return Dict({**{f"{name}":data},
                     **payload})


class OutputProcess(Process):
    def __call__(self, *, payload:Dict, outputs:dict, **kwargs) -> Dict:
        """Output process to save data to different destinations.

        Parameters
        ----------
        payload : Dict
            The data payload to be saved.
        outputs : dict
            A dictionary where keys are names in the payload and values are
            dictionaries with 'package', 'method', and 'dest' keys to specify
            where and how to save each item.
        """

        for name, output_path in outputs.items():
            if isinstance(output_path, list):
                payload_item = payload.get(name)
                if payload_item is None:
                        continue
                for i, output_spec in enumerate(output_path):
                    if isinstance(payload_item[i], DataFrame):
                        payload_item[i].to_csv(output_spec)
                        continue

                    raise IncompatibleArgsException(
                        f"OutputProcess cannot handle the type of '{type(payload_item)}' for '{name}'"
                    )

            else:
                payload_item = payload.get(name)
                if payload_item is None:
                    continue

                if isinstance(payload_item, DataFrame):
                    payload_item.to_csv(output_path)
                    continue

                raise IncompatibleArgsException(
                    f"OutputProcess cannot handle the type of '{type(payload_item)}' for '{name}'"
                )


        return payload


def build_pipeline_from_yaml(
    yaml_path: str | Path,
    pipeline_name: str,
    output_dir: Optional[str | Path] = None,
    input_sources: Optional[dict[str, str]] = None
) -> DFPipeline:
    """Build a DFPipeline from a YAML configuration file.

    Parameters
    ----------
    yaml_path : str | Path
        Path to the YAML configuration file
    pipeline_name : str
        Name of the pipeline to build from the YAML file
    output_dir : str | Path, optional
        Output directory to prepend to all output file paths from YAML.
        If None, uses paths as specified in YAML.
    input_sources : dict[str, str], optional
        Dictionary mapping input names to source paths. This overrides the 'src'
        field in the YAML for specified inputs. Allows reusing the same pipeline
        configuration with different input files.
        Example: {'raw_data': 'file1.csv', 'meta_data': 'metadata1.csv'}

    Returns
    -------
    DFPipeline
        A configured DFPipeline ready to execute

    Raises
    ------
    ValueError
        If the YAML structure is invalid or pipeline_name is not found
    FileNotFoundError
        If the YAML file does not exist

    Examples
    --------
    # Use default sources from YAML
    pipeline = build_pipeline_from_yaml('config.yaml', 'pipeline_1')

    # Override input sources
    pipeline = build_pipeline_from_yaml(
        'config.yaml',
        'pipeline_1',
        input_sources={'raw_data': 'data/file1.csv'}
    )

    # Process multiple files with same pipeline
    for file in data_files:
        pipeline = build_pipeline_from_yaml(
            'config.yaml',
            'pipeline_1',
            input_sources={'raw_data': file},
            output_dir=f'results/{file.stem}'
        )
        result = pipeline()
    """
    # Load YAML file
    yaml_path = Path(yaml_path)
    if not yaml_path.exists():
        raise FileNotFoundError(f"YAML file not found: {yaml_path}")

    # Convert output_dir to Path if provided
    if output_dir is not None:
        output_dir = Path(output_dir)

    with open(yaml_path, 'r') as f:
        config = yaml.safe_load(f)

    # Validate top-level structure
    if 'pipelines' not in config:
        raise ValueError("YAML must contain 'pipelines' key")

    pipelines = config['pipelines']

    # Find the requested pipeline
    pipeline_config = None
    for pipeline_dict in pipelines:
        if pipeline_name in pipeline_dict:
            pipeline_config = pipeline_dict[pipeline_name]
            break

    if pipeline_config is None:
        available = [list(p.keys())[0] for p in pipelines]
        raise ValueError(
            f"Pipeline '{pipeline_name}' not found. Available pipelines: {available}"
        )

    # Validate required sections
    if 'Inputs' not in pipeline_config:
        raise ValueError(f"Pipeline '{pipeline_name}' must contain 'Inputs' section")
    if 'Processes' not in pipeline_config:
        raise ValueError(f"Pipeline '{pipeline_name}' must contain 'Processes' section")
    if 'Outputs' not in pipeline_config:
        raise ValueError(f"Pipeline '{pipeline_name}' must contain 'Outputs' section")

    processes = []

    # Build InputProcess for each input
    inputs_config = pipeline_config['Inputs']
    for input_dict in inputs_config:
        for input_name, input_spec in input_dict.items():
            # Convert list of dicts to a single dict
            input_params = {}
            for item in input_spec:
                input_params.update(item)

            # Override src with input_sources if provided
            if input_sources and input_name in input_sources:
                input_params['src'] = input_sources[input_name]

            # Validate required fields
            if 'src' not in input_params:
                raise ValueError(
                    f"Input '{input_name}' must have 'src' field in YAML or provided via input_sources parameter"
                )
            if 'package' not in input_params:
                raise ValueError(f"Input '{input_name}' must have 'package' field")
            if 'method' not in input_params:
                raise ValueError(f"Input '{input_name}' must have 'method' field")

            # Create a ProcessLogic wrapper for InputProcess
            def make_input_process(name, params):
                def input_logic(**kwargs):
                    input_proc = InputProcess()
                    return input_proc(
                        name=name,
                        src=params['src'],
                        package=params['package'],
                        method=params['method'],
                        **{k: v for k, v in params.items() if k not in ['src', 'package', 'method']}
                    )
                return input_logic

            processes.append(ProcessLogic(make_input_process(input_name, input_params)))

    # Build DFProcess for each process
    processes_config = pipeline_config['Processes']
    for process_dict in processes_config:
        for process_name, process_spec in process_dict.items():
            # Validate required fields
            if 'package' not in process_spec:
                raise ValueError(f"Process '{process_name}' must have 'package' field")
            if 'method' not in process_spec:
                raise ValueError(f"Process '{process_name}' must have 'method' field")
            if 'parameters' not in process_spec:
                raise ValueError(f"Process '{process_name}' must have 'parameters' field")

            # Create a ProcessLogic wrapper for DFProcess
            def make_df_process(name, spec):
                def df_logic(**kwargs):
                    df_proc = DFProcess()
                    return df_proc(
                        payload=Dict(**kwargs),
                        name=name,
                        package=spec['package'],
                        method=spec['method'],
                        parameters=spec['parameters']
                    )
                return df_logic

            processes.append(ProcessLogic(make_df_process(process_name, process_spec)))

    # Build OutputProcess
    outputs_config = pipeline_config['Outputs']
    outputs_dict = {}
    for output_dict in outputs_config:
        outputs_dict.update(output_dict)

    # Combine output_dir with output file paths if output_dir is provided
    if output_dir is not None:
        processed_outputs = {}
        for key, value in outputs_dict.items():
            if isinstance(value, list):
                # Handle list of output paths
                processed_outputs[key] = [str(output_dir / v) for v in value]
            else:
                # Handle single output path
                processed_outputs[key] = str(output_dir / value)
        outputs_dict = processed_outputs

    def output_logic(**kwargs):
        output_proc = OutputProcess()
        return output_proc(payload=Dict(**kwargs), outputs=outputs_dict)

    processes.append(ProcessLogic(output_logic))

    # Create and return the DFPipeline
    return DFPipeline(processes)
