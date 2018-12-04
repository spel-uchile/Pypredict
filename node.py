
class Node:

    def __init__(self, name="", lat=0, lng=0, alt=0, freq=437225000):
        self.name = name
        self.lat = lat
        self.lng = lng
        self.alt = alt
        self.freq = freq

    def __call__(self):
        return self

    def __str__(self):
        return self.name

    def setLat(self, lat):
        self.lat = lat

    def setLng(self, lng):
        self.lng = lng

    def setAlt(self, alt):
        self.alt = alt

    def setFreq(self, freq):
        self.freq = freq

    def getLat(self):
        return self.lat

    def getLng(self):
        return self.lng

    def getAlt(self):
        return self.alt

    def getFreq(self):
        return self.freq
