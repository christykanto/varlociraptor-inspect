from varlociraptor_inspect.views.main import main_view

await main_view()  # noqa: F704, PLE1142  # pyrefly: ignore[invalid-syntax] — stlite compiles with top-level await support
