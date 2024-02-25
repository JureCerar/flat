if __name__.endswith(".fullinit"):
    from . import init
   
    init()

else:
    from pathlib import Path
    from importlib import util
    import sys

    try:
        _script_path = __script__  # type: ignore[name-defined]
    except NameError:
        raise ImportError("invalid invocation of `flat.fullinit`") from None

    path = Path(_script_path).parent / "__init__.py"
    spec = util.spec_from_file_location("flat", path)
    assert spec is not None
    assert spec.loader is not None
    module = util.module_from_spec(spec)
    assert module is not None
    sys.modules["flat"] = module
    spec.loader.exec_module(module)

    import flat.fullinit
