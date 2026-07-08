# log_config.py
import logging

def configure_logger(name, level=logging.INFO) -> logging.Logger:
    logger = logging.getLogger(name)

    logging.basicConfig(
        level=level,
        format='%(asctime)s %(levelname)s [%(name)s] %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )

    return logger