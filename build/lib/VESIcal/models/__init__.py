from VESIcal.models import shishkina
from VESIcal.models import dixon
from VESIcal.models import iaconomarziano
from VESIcal.models import liu
from VESIcal.models import moore
from VESIcal.models import allison

default_models = {'ShishkinaIdealMixing':     shishkina.mixed,
                  'Dixon':                    dixon.mixed,
                  'IaconoMarziano':           iaconomarziano.mixed,
                  'Liu':                      liu.mixed,
                  'ShishkinaCarbon':          shishkina.carbon(),
                  'ShishkinaWater':           shishkina.water(),
                  'DixonCarbon':              dixon.carbon(),
                  'DixonWater':               dixon.water(),
                  'IaconoMarzianoCarbon':     iaconomarziano.carbon(),
                  'IaconoMarzianoWater':      iaconomarziano.water(),
                  'AllisonCarbon':            allison.vesuvius,
                  'AllisonCarbon_sunset':     allison.sunset,
                  'AllisonCarbon_sfvf':       allison.sfvf,
                  'AllisonCarbon_erebus':     allison.erebus,
                  'AllisonCarbon_vesuvius':   allison.vesuvius,
                  'AllisonCarbon_etna':       allison.etna,
                  'AllisonCarbon_stromboli':  allison.stromboli,
                  'MooreWater':               moore.water(),
                  'LiuWater':                 liu.water(),
                  'LiuCarbon':                liu.carbon()
                  }


def get_models(model='all'):
    """
    Returns model objects as a list

    Parameters
    ----------
    models:    str
        OPTIONAL. Default value is 'all' in which case all keys in default_models are returned.
        If 'mixed' is passed, only the MixedFluid model names are returned.
    """
    if model == 'all':
        return list(default_models.values())
    if model == 'mixed':
        # MagmaSat not included here as it is treated separately
        return [shishkina.mixed, dixon.mixed, iaconomarziano.mixed, liu.mixed]


def get_model_names(model='all'):
    """
    Returns all available model names as a list of strings.

    Parameters
    ----------
    models:    str
        OPTIONAL. Default value is 'all' in which case all keys in default_models are returned.
        If 'mixed' is passed, only the MixedFluid model names are returned.
    """
    if model == 'all':
        return list(default_models.keys())
    if model == 'mixed':
        # MagmaSat not included here as it is treated separately
        return ['ShishkinaIdealMixing', 'Dixon', 'IaconoMarziano', 'Liu']
