"""Public package exports for shared fab commands."""

from .project import DEFAULT_HOSTS
from .tasks import register_tasks

__all__ = ["DEFAULT_HOSTS", "register_tasks"]
