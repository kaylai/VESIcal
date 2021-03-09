from VESIcal.models import magmasat
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
